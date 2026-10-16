"""The standard metrics in eval/metrics.py, checked against hand-worked values."""

import json
from pathlib import Path

import pytest

from eval import run
from eval.metrics import (citation_checks, context_scores, grounding, hit_at_k, identifiers,
                          list_scores, ndcg_at_k, precision_recall, reciprocal_rank, wilson_ci)

GENES = frozenset({"PPP1CA", "EVL", "NRAS", "TP53", "WAS", "MAP2K1", "PC"})
PPP1CA = {"source_ids": ["depmap:row648", "peptide:row327", "ms:row1"], "gene": "PPP1CA",
          "protein_change": "p.D242N", "peptide": {"peptide_sequence_fwd": "KFLHKHDLNLICRAHQVVEDGYEFFAK"},
          "ms_detection": [{"peptide": "AHQVVEDGYEFFAK", "covers_variant": False}]}
EVL = {"source_ids": ["depmap:row1139", "ms:row3"], "gene": "EVL", "protein_change": "p.D20N",
       "ms_detection": [{"peptide": "ASVMVYDNTSK", "covers_variant": True}]}


def test_ranking_metrics():
    ranked = ["a", "x", "b"]
    assert hit_at_k(ranked, {"b"}, 1) == 0 and hit_at_k(ranked, {"b"}, 5) == 1
    assert reciprocal_rank(ranked, {"b"}) == pytest.approx(1 / 3)
    assert reciprocal_rank(ranked, {"zzz"}) == 0
    assert precision_recall(ranked, {"a", "b", "c", "d"}) == (pytest.approx(2 / 3), 0.5)
    # relevant at ranks 1 and 3: (1 + 1/log2(4)) / (1 + 1/log2(3))
    assert ndcg_at_k(ranked, {"a", "b"}, 10) == pytest.approx(0.9197, abs=1e-4)
    assert ndcg_at_k(["a", "b"], {"a", "b"}, 10) == 1.0


def test_wilson_interval():
    assert wilson_ci(8, 8) == (pytest.approx(0.6756, abs=1e-4), 1.0)
    assert wilson_ci(0, 10) == (0.0, pytest.approx(0.2775, abs=1e-4))
    low, high = wilson_ci(5, 10)
    assert low == pytest.approx(0.2366, abs=1e-4) and high == pytest.approx(0.7634, abs=1e-4)
    assert wilson_ci(0, 0) is None


def test_identifiers():
    found = identifiers("PPP1CA p.Asp242Asn (rs123) has peptide AHQVVEDGYEFFAK [VCV000013901]; "
                        "MAP2K1 was in ClinVar, CLINVAR and MISSENSE are words, PC too", GENES)
    assert found == {("gene", "PPP1CA"), ("gene", "MAP2K1"), ("change", "D242N"), ("rsid", "rs123"),
                     ("peptide", "AHQVVEDGYEFFAK"), ("vcv", "13901")}  # no "was" (WAS), "P2K" or PC


def test_grounding_flags_planted_hallucinations():
    good = "PPP1CA D242N has an MS hit, AHQVVEDGYEFFAK, that does not span the mutation."
    assert grounding(good, "Was PPP1CA D242N detected?", [PPP1CA], GENES)["ungrounded"] == []
    bad = "PPP1CA D242N was detected as AHQVVEDGYEFFAR (rs999); see also EVL [clinvar:VCV42]."
    assert grounding(bad, "Was PPP1CA D242N detected?", [PPP1CA], GENES)["ungrounded"] == [
        "AHQVVEDGYEFFAR", "EVL", "VCV42", "rs999"]
    # a peptide that is part of a longer sequence in the evidence is grounded
    assert grounding("KFLHKHDLNLICR", "", [PPP1CA], GENES)["ungrounded"] == []


def test_citation_precision_and_recall():
    answer = ("PPP1CA D242N has the MS peptide AHQVVEDGYEFFAK [ms:row1]. "
              "EVL D20N was detected as ASVMVYDNTSK [ms:row1]. "
              "NRAS G12C is a hotspot. "
              "KRAS was not found in the data.")
    c = citation_checks(answer, [PPP1CA, EVL], GENES, skip=lambda s: "not found" in s)
    # 2 citations; the second points to PPP1CA's record, which says nothing about EVL
    assert (c["cited"], c["relevant"], c["unsupported_citations"]) == (2, 1, ["ms:row1"])
    # 3 claims (the "not found" sentence is skipped); only the first is fully supported
    assert (c["claims"], c["supported_claims"]) == (3, 1)


def test_citation_recall_needs_every_identifier():
    # rs999 is not in the cited record, even though D242N is
    c = citation_checks("PPP1CA D242N is rs999 [depmap:row648].", [PPP1CA], GENES)
    assert (c["relevant"], c["supported_claims"]) == (1, 0)


def test_list_and_context_scores():
    s = list_scores("EVL, PPP1CA and TP53 are listed.", ["EVL", "TP53", "NRAS"], "genes", GENES)
    assert (s["precision"], s["recall"]) == (pytest.approx(2 / 3), pytest.approx(2 / 3))
    assert s["missing"] == ["NRAS"] and s["extra"] == ["PPP1CA"]
    assert list_scores("Its variants are P2514Rfs*4 and A348T.", ["P2514fs", "A348T"],
                       "changes")["f1"] == 1.0
    assert list_scores("None listed.", ["EVL"], "genes", GENES)["f1"] == 0.0
    ctx = context_scores([PPP1CA, EVL, PPP1CA], ["depmap:row1139", "depmap:row9"])
    assert (ctx["precision"], ctx["recall"], ctx["missed"]) == (0.5, 0.5, ["depmap:row9"])


def test_consistency_across_runs():
    runs = [[{"id": "a", "passed": True}, {"id": "b", "passed": True}, {"id": "c", "passed": True}],
            [{"id": "a", "passed": True}, {"id": "b", "passed": False}, {"id": "c", "error": "429"}]]
    assert run.consistency(runs) == {"runs": 2, "cases": 2, "consistent": 1, "pass_every_run": 1}


def test_calibrate_judge_counts_caught_and_missed_errors():
    items = [{"id": "1", "question": "q", "answer": "a", "evidence": [], "faithful": True, "planted_error": None},
             {"id": "2", "question": "q", "answer": "b", "evidence": [], "faithful": False, "planted_error": "x"},
             {"id": "3", "question": "q", "answer": "c", "evidence": [], "faithful": False, "planted_error": "y"},
             {"id": "4", "question": "q", "answer": "d", "evidence": [], "faithful": True, "planted_error": None}]
    verdicts = {"a": 1.0, "b": 0.5, "c": 1.0, "d": 0.8}  # caught 2, missed 3, flagged faithful 4
    c = run.calibrate_judge(lambda q, a, e: {"faithfulness": verdicts[a]}, items)
    assert (c["accuracy_n"], c["judged"]) == (2, 4)
    assert c["precision"] == 0.5 and c["recall"] == 0.5


def test_calibration_set_is_balanced_and_labelled():
    items = [json.loads(l) for l in open(Path(__file__).parent.parent / "eval" / "judge_calibration.jsonl")]
    assert {i["faithful"] for i in items} == {True, False}
    assert all(i["evidence"] and (i["faithful"] or i["planted_error"]) for i in items)


def test_retrieval_eval_kinds(index):
    questions = [json.loads(l) for l in open(Path(__file__).parent.parent / "eval" / "questions.jsonl")]
    by_style = {}
    for q in questions:
        by_style.setdefault(q["style"], q)
    rows = run.evaluate_retrieval(index, list(by_style.values()))
    kinds = {r["style"]: r for r in rows}
    assert kinds["one_letter"]["hit1"] == 1 and kinds["one_letter"]["rr"] == 1.0
    assert kinds["gene_list"]["recall"] == 1.0 and kinds["gene_list"]["precision"] == 1.0
    assert kinds["absent"]["no_false_match"]


def test_generated_eval_sets_match_the_data(index):
    from eval.make_answers import calibration_items, generated_cases

    eval_dir = Path(__file__).parent.parent / "eval"
    for name, expected in (("answers_generated.jsonl", generated_cases(index)),
                           ("judge_calibration.jsonl", calibration_items(index))):
        on_disk = [json.loads(l) for l in open(eval_dir / name)]
        assert on_disk == json.loads(json.dumps(expected)), f"rerun python -m eval.make_answers ({name})"
