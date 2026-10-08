"""The answer scorecard in eval/: scoring rules, refusal detection, judge parsing, runner."""

import json
from pathlib import Path

import pytest

from eval import run
from eval.scoring import is_refusal, parse_judgement, score_answer
from rag_assistant import config

CASES = [json.loads(line) for line in open(Path(__file__).parent.parent / "eval" / "answers.jsonl")]


def test_answer_cases_are_well_formed():
    assert len({c["id"] for c in CASES}) == len(CASES)
    for c in CASES:
        assert c["category"] in {"fact", "list", "trap", "refusal"}
        assert c.get("must_include") or c.get("expect_refusal")


def test_score_passes_correct_cited_answer():
    case = {"must_include": [["TP53"], ["R306", "Arg306"]], "expected_source": "depmap:row1535"}
    result = {"answer": "rs121913344 is TP53 p.Arg306Ter [depmap:row1535].",
              "cited_sources": ["depmap:row1535"], "unverified_citations": []}
    assert score_answer(case, result)["passed"]


def test_score_reports_each_failure():
    case = {"must_include": [["TP53"], ["R306"]], "must_not_include": ["sotorasib"],
            "expected_source": "depmap:row1535"}
    result = {"answer": "TP53 can be targeted by sotorasib [clinvar:VCV1].",
              "cited_sources": ["clinvar:VCV1"], "unverified_citations": ["clinvar:VCV1"]}
    s = score_answer(case, result)
    assert not s["passed"] and not s["correct"]
    assert s["missing"] == ["R306"] and s["forbidden"] == ["sotorasib"]
    assert s["cited_expected"] is False and s["unverified_citations"] == 1


@pytest.mark.parametrize("text, refused", [
    ("BRAF V600E was not found in the indexed MOLT4 data.", True),
    ("MOLT4 does not have a KRAS mutation in these records.", True),
    ("There is no information on drugs in the retrieved records.", True),
    ("NRAS G12C is a pathogenic hotspot [depmap:row183].", False),
])
def test_is_refusal(text, refused):
    assert is_refusal(text) is refused


def test_refusal_case_needs_refusal():
    case = {"expect_refusal": True}
    assert score_answer(case, {"answer": "Not found in the indexed data."})["correct"]
    assert not score_answer(case, {"answer": "Yes, BRAF V600E is present."})["correct"]


def test_parse_judgement():
    assert parse_judgement('{"supported": 3, "unsupported": 1, "unsupported_claims": ["x"]}') == {
        "faithfulness": 0.75, "unsupported_claims": ["x"]}
    assert parse_judgement('Sure:\n```json\n{"supported": 2, "unsupported": 0}\n```')["faithfulness"] == 1.0
    assert parse_judgement("no json here") is None


def test_run_answers_and_summary():
    answers = {
        "q1": {"answer": "TP53 R306Ter [depmap:row1535]", "cited_sources": ["depmap:row1535"],
               "unverified_citations": [], "model": "m", "evidence": []},
        "q2": {"answer": "Yes it is present.", "cited_sources": [], "unverified_citations": [],
               "model": "m", "evidence": []},
    }
    cases = [
        {"id": "a", "category": "fact", "question": "q1", "must_include": [["TP53"]],
         "expected_source": "depmap:row1535"},
        {"id": "b", "category": "refusal", "question": "q2", "expect_refusal": True},
        {"id": "c", "category": "fact", "question": "boom", "must_include": [["x"]]},
    ]

    def fake_ask(question):
        if question == "boom":
            raise RuntimeError("503 high demand")
        return answers[question]

    rows = run.run_answers(cases, fake_ask, judge_fn=lambda q, a, e: {"faithfulness": 0.5})
    summary = run.summarize(rows)
    assert [r["passed"] for r in rows] == [True, False, False]
    assert summary["errors"] == 1 and summary["passed"] == pytest.approx(1 / 3)
    assert summary["by_category"] == {"fact": 0.5, "refusal": 0.0}
    assert summary["faithfulness"] == 0.5


def test_parse_model():
    assert run.parse_model("gemini-3.6-flash") == ("gemini", "gemini-3.6-flash")
    assert run.parse_model("openai/gpt-4o-mini") == ("openai", "gpt-4o-mini")


def test_only_model_restores_settings(monkeypatch):
    monkeypatch.setattr(config, "GEMINI_CHAT_MODEL", "main")
    monkeypatch.setattr(config, "GEMINI_FALLBACK_MODELS", ["backup"])
    with run.only_model("gemini", "tested"):
        assert config.GEMINI_CHAT_MODEL == "tested" and config.GEMINI_FALLBACK_MODELS == []
    assert config.GEMINI_CHAT_MODEL == "main" and config.GEMINI_FALLBACK_MODELS == ["backup"]


def test_judge_uses_configured_model_not_tested_one(monkeypatch):
    from rag_assistant import llm

    monkeypatch.setattr(config, "OPENAI_API_KEY", None)
    monkeypatch.setattr(config, "GEMINI_API_KEY", "gm-test")
    monkeypatch.setattr(config, "GEMINI_CHAT_MODEL", "judge-model")
    monkeypatch.setattr(config, "GEMINI_FALLBACK_MODELS", [])
    used = []

    class Client:
        class chat:
            class completions:
                @staticmethod
                def create(model, **_):
                    used.append(model)
                    msg = type("M", (), {"content": '{"supported": 1, "unsupported": 0}'})
                    return type("R", (), {"choices": [type("C", (), {"message": msg})]})

    monkeypatch.setattr(llm, "_client", lambda key, base_url: Client)
    judge = run.make_judge()
    with run.only_model("gemini", "tested-model"):
        assert judge("q", "a", [])["faithfulness"] == 1.0
        assert config.GEMINI_CHAT_MODEL == "tested-model"  # restored after judging
    assert used == ["judge-model"]
