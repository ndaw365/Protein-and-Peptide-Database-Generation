"""Command-line entry point: python -m rag_assistant <command> ..."""

import argparse
import json
import sys

from . import config


def _print_records(records):
    if not records:
        print("No matching variants in the index.")
    for r in records:
        print(f"\n== {r['gene']} {r['protein_change'] or ''} ({r['variant_info']})  "
              f"[{', '.join(r['source_ids'])}]")
        print(json.dumps({k: v for k, v in r.items() if k not in ("card", "source_ids")
                          and v not in (None, [], {})}, indent=2))


def main(argv=None):
    parser = argparse.ArgumentParser(prog="rag_assistant", description=__doc__)
    parser.add_argument("--db", default=str(config.INDEX_DB), help="index path (default: %(default)s)")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("ingest", help="build the variant index")
    p.add_argument("--clinvar", default=config.CLINVAR_URL,
                   help="ClinVar variant_summary.txt(.gz) path or URL (default: NCBI FTP)")
    p.add_argument("--no-clinvar", action="store_true", help="skip ClinVar")
    p.add_argument("--detections", default=str(config.DETECTIONS_TSV),
                   help="FragPipe combined_peptide.tsv from the variant peptide search")
    p.add_argument("--gene-effect", help="DepMap CRISPRGeneEffect.csv (optional)")
    p.add_argument("--depmap-model", help="MOLT4's DepMap model ID (ACH-...) for --gene-effect")
    p.add_argument("--embed", action="store_true", help="embed variant cards for vector search")

    p = sub.add_parser("embed", help="embed variant cards; resumes where a previous run stopped")
    p.add_argument("--provider", choices=["openai", "gemini"],
                   help="embedding provider (default: the one already used, else OpenAI then Gemini)")

    p = sub.add_parser("lookup", help="exact variant lookup, no LLM")
    p.add_argument("gene", nargs="?")
    p.add_argument("protein_change", nargs="?")
    p.add_argument("--rsid")
    p.add_argument("--uniprot")

    p = sub.add_parser("search", help="keyword/vector search, no LLM")
    p.add_argument("query")
    p.add_argument("-k", type=int, default=5)

    p = sub.add_parser("ask", help="answer a question with an LLM grounded on the index")
    p.add_argument("question")
    p.add_argument("--provider", choices=["openai", "gemini"],
                   help="use only this provider (default: OpenAI, falling back to Gemini)")
    p.add_argument("--no-llm", action="store_true", help="only show what would be retrieved")
    p.add_argument("--json", action="store_true", help="print the full result as JSON")

    args = parser.parse_args(argv)

    if args.command == "ingest":
        from .ingest import build_index

        stats = build_index(
            db_path=args.db, clinvar=None if args.no_clinvar else args.clinvar,
            gene_effect=args.gene_effect, depmap_model=args.depmap_model, embed=args.embed,
            detections=args.detections)
        print(f"Index written to {args.db}")
        for key, value in stats.items():
            print(f"  {key}: {value}")
        return 1 if "embedding_error" in stats else 0

    if args.command == "embed":
        from .ingest import embed_index
        from .llm import NoProviderError, error_hint

        try:
            done_now, total, cards, error = embed_index(args.db, provider=args.provider)
        except NoProviderError as exc:
            print(exc, file=sys.stderr)
            return 2
        print(f"Embedded {done_now} cards this run; {total}/{cards} cards have embeddings.")
        if error:
            print(f"Stopped early: {error}", file=sys.stderr)
            print(error_hint(error) or "Run the same command again to resume.", file=sys.stderr)
            return 1
        return 0

    from .retrieve import Index

    try:
        index = Index(args.db)
    except (FileNotFoundError, RuntimeError) as exc:  # no index yet, or built by an older version
        print(exc, file=sys.stderr)
        return 2
    if args.command == "lookup":
        _print_records(index.lookup(gene=args.gene, protein_change=args.protein_change,
                                    rsid=args.rsid, uniprot=args.uniprot))
    elif args.command == "search":
        _print_records(index.search(args.query, args.k))
    elif args.command == "ask":
        if args.no_llm:
            result = index.retrieve(args.question)
            print(f"Retrieval mode: {result['mode']}  entities: {result['entities']}")
            if result.get("total", 0) > len(result["records"]):
                print(f"Showing {len(result['records'])} of {result['total']} matching variants.")
            _print_records(result["records"])
            return 0
        from .llm import NoProviderError, ProvidersFailed, ask, error_hint

        try:
            result = ask(args.question, index=index, provider=args.provider)
        except NoProviderError as exc:
            print(exc, file=sys.stderr)
            return 2
        except ProvidersFailed as exc:
            print(exc, file=sys.stderr)
            if error_hint(str(exc)):
                print(error_hint(str(exc)), file=sys.stderr)
            return 1
        if args.json:
            print(json.dumps(result, indent=2))
        else:
            print(result["answer"])
            print(f"\n[{result['provider']} {result['model']}, retrieval: {result['retrieval_mode']}]")
            if result["unverified_citations"]:
                print("Warning: cited sources not in the retrieved records: "
                      + ", ".join(result["unverified_citations"]), file=sys.stderr)
    return 0
