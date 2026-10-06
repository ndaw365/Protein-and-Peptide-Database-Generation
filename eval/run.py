"""Measure variant-lookup accuracy and latency of the retriever against a naive baseline.

    python -m rag_assistant ingest ...   # build the index first
    python eval/run.py [--db data/variants.sqlite] [--llm]

The baseline is what you would do without the index: re-read "MOLT4 mutations.csv" and
keep rows whose gene symbol and protein change (as DepMap writes them) both appear in
the question. --llm also runs each question through the LLM and checks that the answer
cites the expected DepMap row.
"""

import argparse
import csv
import json
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rag_assistant import config  # noqa: E402
from rag_assistant.retrieve import Index  # noqa: E402


def baseline_lookup(question):
    with open(config.DEPMAP_MUTATIONS_CSV, newline="", encoding="utf-8") as fh:
        hits = []
        for i, r in enumerate(csv.DictReader(fh), start=1):
            change = r["Protein Change"].removeprefix("p.")
            if r["Gene"] and change and r["Gene"] in question.split() and change in question:
                hits.append(f"depmap:row{i}")
        return hits


def timed(fn, *args):
    start = time.perf_counter()
    result = fn(*args)
    return result, (time.perf_counter() - start) * 1000


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", default=str(config.INDEX_DB))
    parser.add_argument("--questions", default=str(Path(__file__).with_name("questions.jsonl")))
    parser.add_argument("--llm", action="store_true", help="also score LLM answers (needs an API key)")
    args = parser.parse_args()

    questions = [json.loads(line) for line in open(args.questions)]
    index = Index(args.db)
    index.retrieve("warm up")  # build the BM25 index once, outside the timings

    rows, by_style = [], {}
    for q in questions:
        result, ms = timed(index.retrieve, q["question"])
        top = result["records"][0]["source_ids"][0] if result["records"] else None
        base, base_ms = timed(baseline_lookup, q["question"])
        row = {"style": q["style"], "rag_ok": top == q["expected_source"], "rag_ms": ms,
               "base_ok": base[:1] == [q["expected_source"]], "base_ms": base_ms}
        if args.llm:
            from rag_assistant.llm import ask

            answer = ask(q["question"], index=index)
            row["llm_ok"] = q["expected_source"] in answer["cited_sources"]
            row["llm_unverified"] = len(answer["unverified_citations"])
        rows.append(row)
        by_style.setdefault(q["style"], []).append(row)

    def pct(rs, key):
        return f"{100 * sum(r[key] for r in rs) / len(rs):5.1f}%"

    print(f"{'style':16} {'n':>3}  {'indexed top-1':>13}  {'CSV-scan top-1':>14}")
    for style, rs in by_style.items():
        print(f"{style:16} {len(rs):3}  {pct(rs, 'rag_ok'):>13}  {pct(rs, 'base_ok'):>14}")
    print(f"{'all':16} {len(rows):3}  {pct(rows, 'rag_ok'):>13}  {pct(rows, 'base_ok'):>14}")
    print(f"\nmedian latency: indexed {statistics.median(r['rag_ms'] for r in rows):.2f} ms, "
          f"CSV scan {statistics.median(r['base_ms'] for r in rows):.2f} ms")
    if args.llm:
        print(f"LLM answers citing the expected row: {pct(rows, 'llm_ok')}; "
              f"unverified citations: {sum(r['llm_unverified'] for r in rows)}")


if __name__ == "__main__":
    main()
