"""Score LLM answers against the curated cases in eval/answers.jsonl.

Each case can set:
    must_include      list of groups; every group must match, a group is a list of
                      alternatives (case-insensitive substrings)
    must_not_include  substrings that must not appear
    expected_source   a source ID the answer should cite
    expect_refusal    true when the data has no answer and the assistant should say so
"""

import json
import re

_REFUSAL_RE = re.compile(
    r"not found|no (record|data|information|variant|mutation|evidence|such)|not (in|present|"
    r"listed|available|included|contain)|does not (have|contain|include|appear)|"
    r"doesn't (have|contain|include|appear)|isn't (in|listed|present)|no .{0,40} (in|for) (the )?"
    r"(indexed|retrieved|molt4)|cannot (find|answer)|unable to find",
    re.IGNORECASE)


def is_refusal(answer):
    """True when the answer says the information is not in the data."""
    return bool(_REFUSAL_RE.search(answer or ""))


def score_answer(case, result):
    """Checks for one answer. result is the dict returned by rag_assistant.llm.ask."""
    answer = result.get("answer", "")
    low = answer.lower()
    missing = [group for group in case.get("must_include", [])
               if not any(alt.lower() in low for alt in group)]
    forbidden = [s for s in case.get("must_not_include", []) if s.lower() in low]
    expect_refusal = case.get("expect_refusal", False)
    refused = is_refusal(answer)
    expected_source = case.get("expected_source")
    cited_expected = expected_source is None or expected_source in result.get("cited_sources", [])
    unverified = len(result.get("unverified_citations", []))
    correct = not missing and not forbidden and (refused if expect_refusal else True)
    return {
        "correct": correct,
        "missing": [g[0] for g in missing],
        "forbidden": forbidden,
        "refusal_ok": refused == expect_refusal if expect_refusal else None,
        "cited_expected": cited_expected,
        "unverified_citations": unverified,
        "passed": correct and cited_expected and unverified == 0,
    }


JUDGE_PROMPT = """You check whether an answer is supported by evidence records.
List every factual claim in the ANSWER. A claim is supported only if the EVIDENCE states
it. Statements that something was not found count as supported when the evidence does
not contain it. Reply with JSON only:
{"supported": <int>, "unsupported": <int>, "unsupported_claims": ["..."]}"""


def judge_messages(question, answer, evidence):
    return [
        {"role": "system", "content": JUDGE_PROMPT},
        {"role": "user", "content": f"QUESTION: {question}\n\nANSWER:\n{answer}\n\n"
                                    f"EVIDENCE:\n{json.dumps(evidence)[:60000]}"},
    ]


def parse_judgement(text):
    """Faithfulness in [0, 1] plus the unsupported claims, or None if unparseable."""
    match = re.search(r"\{.*\}", text or "", re.DOTALL)
    if not match:
        return None
    try:
        data = json.loads(match.group(0))
        supported, unsupported = int(data["supported"]), int(data["unsupported"])
    except (ValueError, KeyError, TypeError):
        return None
    total = supported + unsupported
    return {"faithfulness": 1.0 if total == 0 else supported / total,
            "unsupported_claims": data.get("unsupported_claims", [])}
