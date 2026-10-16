"""Evaluate the variant assistant.

    python eval/run.py                      # retrieval only: no key needed
    python eval/run.py --llm                # + answer scorecard on eval/answers.jsonl
    python eval/run.py --llm --judge        # + faithfulness, judged by a second LLM call
    python eval/run.py --llm --models gemini-3.6-flash,gemini-3.5-flash   # compare models
    python eval/run.py --llm --repeats 3    # + consistency across runs
    python eval/run.py --calibrate-judge    # how well the judge catches planted errors

Retrieval uses eval/questions.jsonl (from eval/make_questions.py): hit rate@1/@5 and MRR
on single-variant lookups (with a naive CSV scan as baseline), precision/recall/nDCG on
gene lists, and false exact matches on variants MOLT4 lacks. Answers use the hand-written
cases in eval/answers.jsonl (or --answers eval/answers_generated.jsonl, from
eval/make_answers.py); eval/scoring.py and eval/metrics.py define every check. Rates show
95% Wilson confidence intervals. LLM runs save every answer to eval/results/ (git-ignored).
"""

import argparse
import contextlib
import csv
import json
import re
import statistics
import sys
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from eval.metrics import (hit_at_k, ndcg_at_k, precision_recall,  # noqa: E402
                          reciprocal_rank, wilson_ci)
from eval.scoring import judge_messages, parse_judgement, score_answer  # noqa: E402
from rag_assistant import config  # noqa: E402
from rag_assistant.retrieve import Index  # noqa: E402

EVAL_DIR = Path(__file__).resolve().parent


# ---------------------------------------------------------------- retrieval
def baseline_lookup(question):
    """What you'd do without the index: scan the CSV for gene + change as DepMap writes them."""
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


def evaluate_retrieval(index, questions):
    """One row per question with the standard retrieval metrics for its kind."""
    index.retrieve("warm up")  # build the keyword index once, outside the timings
    rows = []
    for q in questions:
        result, ms = timed(index.retrieve, q["question"])
        ranked = [r["source_ids"][0] for r in result["records"]]
        relevant = set(q["relevant"])
        row = {"style": q["style"], "ms": ms}
        if q["style"] == "absent":
            row["kind"] = "absent"
            row["no_false_match"] = result["mode"] != "exact"
        elif q["style"] == "gene_list":
            row["kind"] = "multi"
            row["precision"], row["recall"] = precision_recall(ranked, relevant)
            row["ndcg"] = ndcg_at_k(ranked, relevant, 10)
        else:
            base, base_ms = timed(baseline_lookup, q["question"])
            row.update(kind="single", hit1=hit_at_k(ranked, relevant, 1),
                       hit5=hit_at_k(ranked, relevant, 5), rr=reciprocal_rank(ranked, relevant),
                       base_hit1=hit_at_k(base, relevant, 1), base_ms=base_ms)
        rows.append(row)
    return rows


def rate(successes, n):
    """A rate with its 95% confidence interval, e.g. "95% [83-99]"."""
    if not n:
        return "  n/a"
    low, high = wilson_ci(successes, n)
    return f"{100 * successes / n:3.0f}% [{100 * low:.0f}-{100 * high:.0f}]"


def print_retrieval(rows):
    single = [r for r in rows if r["kind"] == "single"]
    multi = [r for r in rows if r["kind"] == "multi"]
    absent = [r for r in rows if r["kind"] == "absent"]
    print(f"RETRIEVAL  (eval/questions.jsonl, {len(rows)} questions; 95% CIs in brackets)")
    if single:
        print(f"  Single-answer lookups: {'n':>3}  {'hit@1':<14} {'hit@5':<14} {'MRR':>5}  CSV-scan hit@1")
        groups = {}
        for r in single:
            groups.setdefault(r["style"], []).append(r)
        for style, rs in [*groups.items(), ("all", single)]:
            n = len(rs)
            print(f"    {style:20} {n:3}  {rate(sum(r['hit1'] for r in rs), n):<14} "
                  f"{rate(sum(r['hit5'] for r in rs), n):<14} "
                  f"{statistics.mean(r['rr'] for r in rs):5.2f}  {rate(sum(r['base_hit1'] for r in rs), n)}")
    if multi:
        print(f"  Gene lists (n={len(multi)}): precision {statistics.mean(r['precision'] for r in multi):.2f}"
              f"   recall {statistics.mean(r['recall'] for r in multi):.2f}"
              f"   nDCG@10 {statistics.mean(r['ndcg'] for r in multi):.2f}")
    if absent:
        print(f"  Absent variants (n={len(absent)}): no false exact match "
              f"{rate(sum(r['no_false_match'] for r in absent), len(absent))}")
    line = f"  median latency: index {statistics.median(r['ms'] for r in rows):.2f} ms"
    if single:
        line += f", CSV scan {statistics.median(r['base_ms'] for r in single):.2f} ms"
    print(line + "\n")


# ---------------------------------------------------------------- answers
def parse_model(spec):
    """"gemini-3.6-flash", "gemini/gemini-3.6-flash" or "openai/gpt-4o-mini" -> (provider, model)."""
    provider, _, model = spec.partition("/")
    return (provider, model) if model else ("gemini", provider)


@contextlib.contextmanager
def only_model(provider, model):
    """Use exactly one model, with no fallbacks, so the score belongs to that model."""
    names = {"gemini": ("GEMINI_CHAT_MODEL", "GEMINI_FALLBACK_MODELS"),
             "openai": ("OPENAI_CHAT_MODEL", "OPENAI_FALLBACK_MODELS")}[provider]
    saved = [getattr(config, n) for n in names]
    setattr(config, names[0], model)
    setattr(config, names[1], [])
    try:
        yield
    finally:
        for n, v in zip(names, saved):
            setattr(config, n, v)


def short_error(exc):
    """One line per model tried, e.g. "gemini/gemini-3.6-flash: RateLimitError 429"."""
    lines = [l.strip() for l in str(exc).splitlines() if l.strip()]
    tried = [re.match(r"(\S+/\S+): (\w+)(?:: Error code: (\d+))?", l) for l in lines]
    parts = [" ".join(filter(None, m.groups())) for m in tried if m]
    return "; ".join(parts) if parts else f"{type(exc).__name__}: {str(exc)[:200]}"


def is_rate_limit(exc):
    return "RateLimitError" in str(exc) or "429" in str(exc)


def select_cases(cases, only):
    """Cases whose category or id is in the comma-separated `only`; all cases if empty."""
    if not only:
        return cases
    wanted = {name.strip() for name in only.split(",") if name.strip()}
    known = {c["category"] for c in cases} | {c["id"] for c in cases}
    unknown = sorted(wanted - known)
    if unknown:
        sys.exit(f"Unknown --only name(s): {', '.join(unknown)}. Use a category "
                 f"({', '.join(sorted({c['category'] for c in cases}))}) or a case ID "
                 f"from eval/answers.jsonl.")
    return [c for c in cases if c["category"] in wanted or c["id"] in wanted]


def run_answers(cases, ask_fn, judge_fn=None, delay=0.0, retry_wait=60.0, log=print, sleep=time.sleep,
                genes=frozenset()):
    """Ask every case and print one progress line each.

    delay waits between cases (free tiers allow only a few requests per minute). When a
    case hits a rate limit, it waits retry_wait seconds and tries once more. genes (the
    index's gene symbols) lets the hallucination checks recognise gene names.
    """
    rows = []
    for n, case in enumerate(cases, start=1):
        if n > 1 and delay:
            sleep(delay)
        row = {"id": case["id"], "category": case["category"], "question": case["question"]}
        start = time.perf_counter()
        result = None
        for attempt in (1, 2):
            try:
                result = ask_fn(case["question"])
                break
            except Exception as exc:  # provider errors are scored as failures, not crashes
                if attempt == 1 and retry_wait and is_rate_limit(exc):
                    log(f"  [{n}/{len(cases)}] {case['id']}: rate limited, waiting {retry_wait:.0f} s")
                    sleep(retry_wait)
                    continue
                row.update(error=short_error(exc), passed=False, correct=False)
        row["seconds"] = time.perf_counter() - start
        if result is not None:
            row["answer"] = result["answer"]
            row["model"] = result.get("model")
            row.update(score_answer(case, result, genes))
            if judge_fn:
                row["judgement"] = judge_fn(case["question"], result["answer"],
                                            result.get("evidence", []))
        status = "ERROR" if "error" in row else "PASS" if row["passed"] else "FAIL"
        log(f"  [{n}/{len(cases)}] {case['id']:<16} {status:<5} {row['seconds']:5.1f} s")
        rows.append(row)
    return rows


def summarize(rows):
    """Scorecard numbers for one model. Counts are kept so rates can carry CIs."""
    def share(items, key):
        items = [r for r in items if r.get(key) is not None]
        return None if not items else sum(bool(r[key]) for r in items) / len(items)

    def mean(values):
        values = [v for v in values if v is not None]
        return statistics.mean(values) if values else None

    answered = [r for r in rows if "error" not in r]
    by_cat = {}
    for r in rows:
        by_cat.setdefault(r["category"], []).append(r)
    lists = [r["list"] for r in answered if r.get("list")]
    contexts = [r["context"] for r in answered if r.get("context")]
    grounded = [r["grounding"] for r in answered if r.get("grounding")]
    cites = [r["citations"] for r in answered if r.get("citations")]
    refusal_cases = [r for r in answered if r.get("refusal_ok") is not None]
    answerable = [r for r in answered if r.get("false_refusal") is not None]
    faith = [r["judgement"]["faithfulness"] for r in answered if r.get("judgement")]
    return {
        "cases": len(rows),
        "errors": len(rows) - len(answered),
        "passed": share(rows, "passed"),
        "passed_n": sum(bool(r.get("passed")) for r in rows),
        "passed_of_answered": share(answered, "passed"),
        "correct": share(rows, "correct"),
        "correct_n": sum(bool(r.get("correct")) for r in rows),
        "by_category": {c: share(rs, "correct") for c, rs in sorted(by_cat.items())},
        "key_fact_recall": mean(r.get("key_fact_recall") for r in answered),
        "list_precision": mean(x["precision"] for x in lists),
        "list_recall": mean(x["recall"] for x in lists),
        "list_f1": mean(x["f1"] for x in lists),
        "list_cases": len(lists),
        "refusals_ok": sum(r["refusal_ok"] for r in refusal_cases),
        "refusal_cases": len(refusal_cases),
        "false_refusals": sum(r["false_refusal"] for r in answerable),
        "answerable_cases": len(answerable),
        "context_recall": mean(x["recall"] for x in contexts),
        "context_precision": mean(x["precision"] for x in contexts),
        "identifiers": sum(g["identifiers"] for g in grounded),
        "ungrounded": sum(len(g["ungrounded"]) for g in grounded),
        "answers_with_ungrounded": sum(bool(g["ungrounded"]) for g in grounded),
        "answers_checked": len(grounded),
        "citations": sum(c["cited"] for c in cites),
        "citations_relevant": sum(c["relevant"] for c in cites),
        "claims": sum(c["claims"] for c in cites),
        "claims_supported": sum(c["supported_claims"] for c in cites),
        "cited_expected": share(answered, "cited_expected"),
        "unverified_citations": sum(r.get("unverified_citations", 0) for r in answered),
        "faithfulness": statistics.mean(faith) if faith else None,
        "judged": len(faith),
        "median_seconds": statistics.median(r["seconds"] for r in rows) if rows else None,
    }


def consistency(runs):
    """Across repeated runs of the same cases: share with the same outcome every time,
    and share that passed every time (pass^N). Errored attempts are left out."""
    outcomes = {}
    for rows in runs:
        for r in rows:
            if "error" not in r:
                outcomes.setdefault(r["id"], []).append(bool(r["passed"]))
    seen = [o for o in outcomes.values() if len(o) == len(runs)]
    return {"runs": len(runs), "cases": len(seen),
            "consistent": sum(len(set(o)) == 1 for o in seen),
            "pass_every_run": sum(all(o) for o in seen)}


def _pct(x):
    return "  n/a" if x is None else f"{100 * x:3.0f}%"


def print_scorecard(model, summary, rows, consistent=None):
    s = summary
    print(f"ANSWERS  {model}  ({s['cases']} cases; 95% CIs in brackets)")
    print(f"  Overall       passed {rate(s['passed_n'], s['cases'])}   of answered "
          f"{_pct(s['passed_of_answered'])}   errors {s['errors']}   "
          f"median {s['median_seconds']:.1f} s")
    print(f"  Correctness   correct {rate(s['correct_n'], s['cases'])}   key-fact recall "
          f"{_pct(s['key_fact_recall'])}   "
          + "  ".join(f"{c} {_pct(v)}" for c, v in s["by_category"].items()))
    if s["list_cases"]:
        print(f"  Lists         precision {_pct(s['list_precision'])}   recall "
              f"{_pct(s['list_recall'])}   F1 {_pct(s['list_f1'])}   (n={s['list_cases']})")
    print(f"  Abstention    says 'not found' when it should {rate(s['refusals_ok'], s['refusal_cases'])}"
          f"   wrongly says 'not found' {rate(s['false_refusals'], s['answerable_cases'])}")
    print(f"  Context       recall {_pct(s['context_recall'])}   precision "
          f"{_pct(s['context_precision'])}   (relevant records the model saw)")
    print(f"  Grounding     identifiers found in evidence "
          f"{rate(s['identifiers'] - s['ungrounded'], s['identifiers'])}   answers with an "
          f"ungrounded identifier {s['answers_with_ungrounded']}/{s['answers_checked']}")
    print(f"  Citations     precision {rate(s['citations_relevant'], s['citations'])}   recall "
          f"{rate(s['claims_supported'], s['claims'])}   expected source cited "
          f"{_pct(s['cited_expected'])}   unverified {s['unverified_citations']}")
    if s["judged"]:
        print(f"  Faithfulness  judge: {_pct(s['faithfulness'])} of claims supported (n={s['judged']})")
    if consistent:
        c = consistent
        print(f"  Consistency   same outcome in all {c['runs']} runs {rate(c['consistent'], c['cases'])}"
              f"   passed every run {rate(c['pass_every_run'], c['cases'])}")
    for r in rows:
        if not r.get("passed"):
            why = r.get("error") or ", ".join(filter(None, [
                r.get("missing") and f"missing {r['missing']}",
                r.get("forbidden") and f"contains {r['forbidden']}",
                r.get("refusal_ok") is False and "did not say 'not found'",
                r.get("cited_expected") is False and "expected source not cited",
                r.get("unverified_citations") and f"{r['unverified_citations']} unverified citation(s)",
                r.get("grounding", {}).get("ungrounded") and f"ungrounded {r['grounding']['ungrounded']}",
            ]))
            print(f"    FAIL {r['id']}: {why}")
    print()


def calibrate_judge(judge, items):
    """How well the judge flags answers with planted errors (positive = unfaithful)."""
    rows = []
    for item in items:
        verdict = judge(item["question"], item["answer"], item["evidence"])
        flagged = None if verdict is None else verdict["faithfulness"] < 1
        rows.append({"id": item["id"], "faithful": item["faithful"], "flagged": flagged,
                     "planted_error": item["planted_error"]})
    judged = [r for r in rows if r["flagged"] is not None]
    tp = sum(r["flagged"] and not r["faithful"] for r in judged)
    fp = sum(r["flagged"] and r["faithful"] for r in judged)
    fn = sum(not r["flagged"] and not r["faithful"] for r in judged)
    return {"items": len(items), "judged": len(judged),
            "accuracy_n": sum(r["flagged"] != r["faithful"] for r in judged),
            "precision": tp / (tp + fp) if tp + fp else None,
            "recall": tp / (tp + fn) if tp + fn else None, "rows": rows}


def print_calibration(c):
    print(f"JUDGE CALIBRATION  (eval/judge_calibration.jsonl, {c['items']} labelled answers)")
    print(f"  accuracy {rate(c['accuracy_n'], c['judged'])}   catches planted errors (recall) "
          f"{_pct(c['recall'])}   flags are right (precision) {_pct(c['precision'])}")
    for r in c["rows"]:
        if r["flagged"] is None:
            print(f"    {r['id']}: judge gave no usable verdict")
        elif r["flagged"] == r["faithful"]:
            what = f"missed: {r['planted_error']}" if not r["faithful"] else "flagged a faithful answer"
            print(f"    {r['id']}: {what}")
    print()


MODEL_SETTINGS = ("OPENAI_CHAT_MODEL", "OPENAI_FALLBACK_MODELS",
                  "GEMINI_CHAT_MODEL", "GEMINI_FALLBACK_MODELS")


def make_judge():
    """Faithfulness judge: a plain chat call with the models configured in .env.

    The judge keeps those models even while --models swaps the model being tested, so no
    model grades its own answers.
    """
    from rag_assistant.llm import _with_fallback

    judge_settings = {n: getattr(config, n) for n in MODEL_SETTINGS}

    def judge(question, answer, evidence):
        def call(client, _name, chat_model, _embed):
            resp = client.chat.completions.create(
                model=chat_model, messages=judge_messages(question, answer, evidence), temperature=0)
            return resp.choices[0].message.content

        tested = {n: getattr(config, n) for n in MODEL_SETTINGS}
        for n, v in judge_settings.items():
            setattr(config, n, v)
        try:
            return parse_judgement(_with_fallback(call))
        except Exception:
            return None
        finally:
            for n, v in tested.items():
                setattr(config, n, v)
    return judge


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default=str(config.INDEX_DB))
    parser.add_argument("--questions", default=str(EVAL_DIR / "questions.jsonl"))
    parser.add_argument("--answers", default=str(EVAL_DIR / "answers.jsonl"))
    parser.add_argument("--llm", action="store_true", help="score LLM answers (needs an API key)")
    parser.add_argument("--models", help="comma-separated models to compare, e.g. "
                                         "gemini-3.6-flash,openai/gpt-4o-mini")
    parser.add_argument("--judge", action="store_true", help="also judge faithfulness with an LLM")
    parser.add_argument("--only", help="comma-separated categories (fact,list,trap,refusal) "
                                       "and/or case IDs to run, e.g. trap,braf-absent")
    parser.add_argument("--delay", type=float, default=0.0,
                        help="seconds to wait between cases, to stay under per-minute quotas")
    parser.add_argument("--retry-wait", type=float, default=60.0,
                        help="seconds to wait before retrying a rate-limited case once (0 = no retry)")
    parser.add_argument("--repeats", type=int, default=1,
                        help="ask every case N times and report consistency (costs N times the quota)")
    parser.add_argument("--calibrate-judge", action="store_true",
                        help="check the faithfulness judge on eval/judge_calibration.jsonl")
    args = parser.parse_args()

    try:
        index = Index(args.db)
    except (FileNotFoundError, RuntimeError) as exc:  # no index yet, or built by an older version
        sys.exit(str(exc))
    print_retrieval(evaluate_retrieval(index, [json.loads(l) for l in open(args.questions)]))
    if not (args.llm or args.calibrate_judge):
        return

    from rag_assistant.llm import NoProviderError, _providers, ask

    try:
        _providers()
    except NoProviderError as exc:
        sys.exit(str(exc))
    report = {}
    if args.calibrate_judge:
        items = [json.loads(l) for l in open(EVAL_DIR / "judge_calibration.jsonl")]
        report["judge_calibration"] = calibrate_judge(make_judge(), items)
        print_calibration(report["judge_calibration"])
    if args.llm:
        cases = select_cases([json.loads(l) for l in open(args.answers)], args.only)
        judge = make_judge() if args.judge else None
        models = [parse_model(m) for m in args.models.split(",")] if args.models else [(None, None)]
        for provider, model in models:
            label = f"{provider}/{model}" if model else "configured model (with fallbacks)"
            runs = []
            for n in range(1, args.repeats + 1):
                run_label = f" (run {n}/{args.repeats})" if args.repeats > 1 else ""
                print(f"ANSWERS  {label}{run_label}: running {len(cases)} cases")
                with only_model(provider, model) if model else contextlib.nullcontext():
                    runs.append(run_answers(cases, lambda q: ask(q, index=index, provider=provider),
                                            judge, delay=args.delay, retry_wait=args.retry_wait,
                                            log=lambda line: print(line, flush=True),
                                            genes=index.genes))
                print()
            summary = summarize(runs[0])
            consistent = consistency(runs) if len(runs) > 1 else None
            print_scorecard(label, summary, runs[0], consistent)
            report[label] = {"summary": summary, "consistency": consistent, "rows": runs[0],
                             "repeat_rows": runs[1:]}

    out = EVAL_DIR / "results" / f"{datetime.now():%Y%m%d-%H%M%S}.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"Full answers saved to {out.relative_to(EVAL_DIR.parent)}")


if __name__ == "__main__":
    main()
