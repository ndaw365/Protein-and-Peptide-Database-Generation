"""Answer questions with an LLM grounded on the variant index.

OpenAI is tried first. If no OPENAI_API_KEY is set, or the OpenAI call fails
(network, auth, quota, server or model errors), the same request is retried on
Gemini through its OpenAI-compatible endpoint, so both share one code path.
"""

import json
import re

from . import config
from .retrieve import Index

SYSTEM_PROMPT = """You are a variant-lookup assistant for the MOLT4 T-ALL cell line.
Answer ONLY from the retrieved records and tool results. Each fact in them carries a
source ID such as depmap:row12, peptide:row5, clinvar:VCV13901, ms:row3 or dropped:row1.

Rules:
- Cite the source ID in square brackets after each fact, e.g. "ClinVar lists it as
  pathogenic [clinvar:VCV13901]".
- If the records do not contain the answer, say it was not found in the indexed
  data. Never fill gaps from general knowledge.
- "Vep Clin Sig" is DepMap's copy of ClinVar; "clinvar" entries come from the ClinVar
  release that was indexed. Say which one you used.
- An MS hit supports a variant only when covers_variant is true; otherwise the
  detected peptide is identical to the wild-type protein.
- If "total" is larger than the number of records, say the list is partial and how many
  exist; use filter_variants or lookup_variant to fetch more when the question needs them.
- Use the tools to look up more variants when the provided records are not enough.
- Be concise."""

TOOLS = [
    {"type": "function", "function": {
        "name": "lookup_variant",
        "description": "Exact lookup of MOLT4 variants by gene, protein change, dbSNP rsID or UniProt accession.",
        "parameters": {"type": "object", "properties": {
            "gene": {"type": "string", "description": "HGNC symbol, e.g. NRAS"},
            "protein_change": {"type": "string", "description": "e.g. G12C, p.Gly12Cys, R306Ter"},
            "rsid": {"type": "string", "description": "e.g. rs121913250"},
            "uniprot": {"type": "string", "description": "e.g. P01111"},
        }},
    }},
    {"type": "function", "function": {
        "name": "filter_variants",
        "description": ("List MOLT4 variants matching structured filters. Annotation values match "
                        "whole terms: 'pathogenic' does not include 'likely_pathogenic'; call "
                        "twice to get both."),
        "parameters": {"type": "object", "properties": {
            "gene": {"type": "string"},
            "variant_type": {"type": "string", "description": "e.g. missense_variant, frameshift_variant, stop_gained"},
            "clinvar_significance": {"type": "string", "description": "e.g. pathogenic, likely_benign, uncertain_significance"},
            "am_class": {"type": "string", "description": "AlphaMissense class: likely_pathogenic, ambiguous, likely_benign"},
            "flag": {"type": "string", "enum": ["Hotspot", "Hess Driver", "Likely LOF", "Oncogene High Impact",
                                                "Tumor Suppressor High Impact"]},
            "in_peptide_database": {"type": "boolean"},
            "limit": {"type": "integer", "default": 20},
        }},
    }},
    {"type": "function", "function": {
        "name": "search_variants",
        "description": "Free-text search over variant summaries (gene names, descriptions, annotations).",
        "parameters": {"type": "object", "properties": {
            "query": {"type": "string"}, "k": {"type": "integer", "default": 5},
        }, "required": ["query"]},
    }},
]


class NoProviderError(RuntimeError):
    pass


def _providers(preferred=None):
    available = {
        "openai": (config.OPENAI_API_KEY, None, config.OPENAI_CHAT_MODEL, config.OPENAI_EMBED_MODEL),
        "gemini": (config.GEMINI_API_KEY, config.GEMINI_BASE_URL, config.GEMINI_CHAT_MODEL,
                   config.GEMINI_EMBED_MODEL),
    }
    order = [preferred] if preferred else ["openai", "gemini"]
    usable = [(name, *available[name]) for name in order if available[name][0]]
    if not usable:
        wanted = preferred or "openai or gemini"
        raise NoProviderError(
            f"No API key for {wanted}. Set OPENAI_API_KEY and/or GEMINI_API_KEY "
            "(see .env.example), or use --no-llm to see the retrieved records only.")
    return usable


def _client(api_key, base_url):
    from openai import OpenAI

    return OpenAI(api_key=api_key, base_url=base_url, max_retries=2, timeout=60)


def _with_fallback(fn, preferred=None, client_factory=_client):
    """Run fn(client, name, chat_model, embed_model) on each provider until one succeeds."""
    import openai

    errors = []
    for name, key, base_url, chat_model, embed_model in _providers(preferred):
        try:
            return fn(client_factory(key, base_url), name, chat_model, embed_model)
        except openai.APIError as exc:  # covers connection, auth, rate-limit and server errors
            errors.append(f"{name}: {type(exc).__name__}: {exc}")
    raise RuntimeError("All LLM providers failed:\n  " + "\n  ".join(errors))


def embed_texts(texts, provider=None, client_factory=_client, batch_size=100):
    """Return (vectors, provider_name)."""
    def run(client, name, _chat, embed_model):
        vectors = []
        for i in range(0, len(texts), batch_size):
            resp = client.embeddings.create(model=embed_model, input=texts[i:i + batch_size])
            vectors += [d.embedding for d in resp.data]
        return vectors, name

    return _with_fallback(run, provider, client_factory)


def _compact(record):
    return {k: v for k, v in record.items() if k != "card" and v not in (None, [], {})}


def run_tool(index, name, args):
    if name == "lookup_variant":
        records = index.lookup(**{k: v for k, v in args.items() if v})
    elif name == "filter_variants":
        records = index.filter_variants(**args)
    elif name == "search_variants":
        records = index.search(args["query"], args.get("k", 5))
    else:
        return {"error": f"unknown tool {name}"}
    return {"records": [_compact(r) for r in records]}


def _source_ids(payload):
    return set(re.findall(r"\b(?:depmap|peptide|clinvar|ms|dropped|gene_effect):[\w-]+",
                          json.dumps(payload)))


def ask(question, index=None, provider=None, client_factory=_client, max_steps=6):
    """Answer a question. Returns the answer, the provider used, and citation checks."""
    index = index or Index()
    context = index.retrieve(question)
    context_payload = {"retrieval_mode": context["mode"], "entities": context["entities"],
                       "records": [_compact(r) for r in context["records"]]}
    if "total" in context:
        context_payload["total"] = context["total"]

    def run(client, name, chat_model, _embed):
        messages = [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": f"Question: {question}\n\nRetrieved records:\n"
                                        f"{json.dumps(context_payload, indent=1)}"},
        ]
        retrieved = _source_ids(context_payload)
        tool_log = []
        for _ in range(max_steps):
            resp = client.chat.completions.create(model=chat_model, messages=messages,
                                                  tools=TOOLS, temperature=0)
            msg = resp.choices[0].message
            if not msg.tool_calls:
                answer = msg.content or ""
                cited = set(re.findall(r"\[((?:depmap|peptide|clinvar|ms|dropped|gene_effect):[\w-]+)\]",
                                       answer))
                return {
                    "answer": answer, "provider": name, "model": chat_model,
                    "retrieval_mode": context["mode"], "tool_calls": tool_log,
                    "cited_sources": sorted(cited),
                    "unverified_citations": sorted(cited - retrieved),
                }
            messages.append({"role": "assistant", "content": msg.content or "", "tool_calls": [
                {"id": c.id, "type": "function",
                 "function": {"name": c.function.name, "arguments": c.function.arguments}}
                for c in msg.tool_calls]})
            for call in msg.tool_calls:
                try:
                    args = json.loads(call.function.arguments or "{}")
                    result = run_tool(index, call.function.name, args)
                except (TypeError, ValueError) as exc:
                    result = {"error": str(exc)}
                retrieved |= _source_ids(result)
                tool_log.append({"tool": call.function.name, "arguments": call.function.arguments})
                messages.append({"role": "tool", "tool_call_id": call.id,
                                 "content": json.dumps(result)})
        return {"answer": "Stopped: too many tool calls without a final answer.", "provider": name,
                "model": chat_model, "retrieval_mode": context["mode"], "tool_calls": tool_log,
                "cited_sources": [], "unverified_citations": []}

    return _with_fallback(run, provider, client_factory)
