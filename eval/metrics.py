"""Standard RAG evaluation metrics, as plain functions with no I/O.

Retrieval: hit rate@k, MRR, precision/recall, nDCG@k.
Answers:   list precision/recall/F1, context precision/recall, and two deterministic
           hallucination checks that need no LLM judge:
             grounding  - every identifier in the answer (protein change, rsID, peptide
                          sequence, ClinVar VCV ID, MOLT4 gene) must appear in the
                          evidence the model was shown, or in the question
             citations  - each cited source must be a record the model saw that is
                          relevant to its sentence (precision), and each sentence's
                          identifiers must all be in the records it cites (recall)
Reliability: Wilson 95% confidence intervals for rates.
"""

import json
import math
import re

from rag_assistant.normalize import find_protein_changes, find_rsids

SOURCE_RE = re.compile(r"\b(?:depmap|peptide|clinvar|ms|dropped|gene_effect):[\w-]+")
_VCV_RE = re.compile(r"\bVCV0*(\d+)\b")
_WORD_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9-]*[A-Za-z0-9]|[A-Za-z0-9]")
# Upper-case runs of amino-acid letters, e.g. ASVMVYDNTSK. Real words that happen to use
# only those letters are skipped.
_PEPTIDE_RE = re.compile(r"\b[ACDEFGHIKLMNPQRSTVWY]{7,}\b")
_NOT_PEPTIDES = {"CLINVAR", "MISSENSE", "ALPHAMISSENSE", "FRAMESHIFT", "DELETERIOUS",
                 "PEPTIDES", "VARIANTS"}


# ---------------------------------------------------------------- retrieval
def hit_at_k(ranked, relevant, k):
    """1 if any relevant item is in the top k, else 0."""
    return int(any(r in relevant for r in ranked[:k]))


def reciprocal_rank(ranked, relevant):
    """1 / rank of the first relevant item (0 if none); its mean is MRR."""
    for rank, r in enumerate(ranked, start=1):
        if r in relevant:
            return 1 / rank
    return 0.0


def precision_recall(returned, relevant):
    """Share of returned items that are relevant, and share of relevant items returned."""
    hits = len(set(returned) & set(relevant))
    precision = hits / len(set(returned)) if returned else 0.0
    recall = hits / len(set(relevant)) if relevant else 0.0
    return precision, recall


def ndcg_at_k(ranked, relevant, k):
    """Normalised discounted cumulative gain with binary relevance."""
    dcg = sum(1 / math.log2(i + 1) for i, r in enumerate(ranked[:k], start=1) if r in relevant)
    ideal = sum(1 / math.log2(i + 1) for i in range(1, min(k, len(relevant)) + 1))
    return dcg / ideal if ideal else 0.0


def wilson_ci(successes, n, z=1.96):
    """95% Wilson score interval for a proportion, as (low, high); None if n is 0."""
    if not n:
        return None
    p = successes / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(0.0, centre - half), min(1.0, centre + half)


def f1(precision, recall):
    return 0.0 if not precision + recall else 2 * precision * recall / (precision + recall)


# ---------------------------------------------------------------- identifiers
def identifiers(text, genes=frozenset()):
    """The checkable identifiers in a text, as a set of (kind, value) pairs."""
    text = text or ""
    found = {("change", c) for c in find_protein_changes(text)}
    found |= {("rsid", r.lower()) for r in find_rsids(text)}
    found |= {("vcv", v) for v in _VCV_RE.findall(text)}
    found |= {("peptide", p) for p in _PEPTIDE_RE.findall(text) if p not in _NOT_PEPTIDES}
    # Case-sensitive, so ordinary words ("was", "set") aren't read as gene symbols, and
    # 1-2 character symbols (PC, HR, TG) are skipped as too likely to be abbreviations
    found |= {("gene", w) for w in _WORD_RE.findall(text) if len(w) > 2 and w in genes}
    return found


class _Vocabulary:
    """What a body of evidence mentions, for fast membership checks."""

    def __init__(self, text):
        self.text = text.upper()
        self.changes = set(find_protein_changes(text))
        self.rsids = {r.lower() for r in find_rsids(text)}
        self.vcvs = set(_VCV_RE.findall(text)) | set(re.findall(r'"variation_id": "0*(\d+)"', text))
        self.words = {w.upper() for w in _WORD_RE.findall(text)}

    def __contains__(self, ident):
        kind, value = ident
        if kind == "change":
            return value in self.changes
        if kind == "rsid":
            return value in self.rsids
        if kind == "vcv":
            return value in self.vcvs
        if kind == "peptide":
            return value in self.text  # may be part of a longer sequence
        return value.upper() in self.words


def _label(ident):
    kind, value = ident
    return f"VCV{value}" if kind == "vcv" else value


def grounding(answer, question, evidence, genes=frozenset()):
    """Identifiers in the answer that are neither in the evidence nor in the question."""
    seen = _Vocabulary(json.dumps(evidence) + "\n" + (question or ""))
    found = identifiers(answer, genes)
    ungrounded = sorted(_label(i) for i in found if i not in seen)
    return {"identifiers": len(found), "ungrounded": ungrounded}


def _sentences(text):
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+|\n+", text or "") if s.strip()]


def citation_checks(answer, evidence, genes=frozenset(), skip=lambda sentence: False):
    """Citation precision and recall, in the sense of the ALCE benchmark.

    A "claim" is a sentence that states at least one identifier; sentences for which
    skip() is true (e.g. "X was not found") are not claims.
    precision: share of citations whose record the model saw and that mentions at least
               one identifier of the citing sentence (the citation is relevant)
    recall:    share of claims whose identifiers all appear in the records they cite
               (an uncited claim counts as unsupported)
    """
    vocab = {}
    for record in evidence:
        text = json.dumps(record)
        for sid in record.get("source_ids", []):
            vocab.setdefault(sid, []).append(_Vocabulary(text))
    cited = relevant = claims = supported_claims = 0
    unsupported_citations, unsupported_claims = [], []
    for sentence in _sentences(answer):
        cites = list(dict.fromkeys(SOURCE_RE.findall(sentence)))
        idents = identifiers(SOURCE_RE.sub("", sentence), genes)
        for sid in cites:
            cited += 1
            ok = sid in vocab and (not idents or any(i in v for v in vocab[sid] for i in idents))
            relevant += ok
            if not ok:
                unsupported_citations.append(sid)
        if idents and not skip(sentence):
            claims += 1
            cited_vocab = [v for sid in cites for v in vocab.get(sid, [])]
            if all(any(i in v for v in cited_vocab) for i in idents):
                supported_claims += 1
            else:
                unsupported_claims.append(sentence[:120])
    return {"cited": cited, "relevant": relevant, "claims": claims,
            "supported_claims": supported_claims, "unsupported_citations": unsupported_citations,
            "unsupported_claims": unsupported_claims}


def list_scores(answer, expected, kind, genes=frozenset()):
    """Precision, recall and F1 of the genes or protein changes an answer lists."""
    if kind == "genes":
        named = {v for k, v in identifiers(answer, genes) if k == "gene"}
    else:
        named = set(find_protein_changes(answer))
    expected = set(expected)
    precision, recall = precision_recall(sorted(named), sorted(expected))
    if not named:
        precision = 0.0
    return {"precision": precision, "recall": recall, "f1": f1(precision, recall),
            "missing": sorted(expected - named), "extra": sorted(named - expected)}


def context_scores(evidence, relevant):
    """Context precision and recall: the records the model saw vs the relevant ones."""
    seen = list(dict.fromkeys(r["source_ids"][0] for r in evidence if r.get("source_ids")))
    precision, recall = precision_recall(seen, relevant)
    return {"precision": precision, "recall": recall, "seen": len(seen),
            "missed": sorted(set(relevant) - set(seen))}
