"""LLM layer tests with fake clients: no network or API keys needed."""

import json
from types import SimpleNamespace

import httpx
import openai
import pytest

from rag_assistant import config, llm


def _message(content=None, tool_calls=None):
    return SimpleNamespace(choices=[SimpleNamespace(message=SimpleNamespace(
        content=content, tool_calls=tool_calls))])


def _tool_call(name, args):
    return SimpleNamespace(id="call_1", function=SimpleNamespace(name=name, arguments=json.dumps(args)))


class FakeClient:
    def __init__(self, replies=None, error=None):
        self.replies, self.error, self.requests = list(replies or []), error, []
        self.chat = SimpleNamespace(completions=SimpleNamespace(create=self._create))

    def _create(self, **kwargs):
        self.requests.append(kwargs)
        if self.error:
            raise self.error
        return self.replies.pop(0)


@pytest.fixture
def both_keys(monkeypatch):
    monkeypatch.setattr(config, "OPENAI_API_KEY", "sk-test")
    monkeypatch.setattr(config, "GEMINI_API_KEY", "gm-test")


def test_falls_back_to_gemini_when_openai_fails(index, both_keys):
    request = httpx.Request("POST", "https://api.openai.com/v1/chat/completions")
    openai_client = FakeClient(error=openai.APIConnectionError(request=request))
    gemini_client = FakeClient(replies=[
        _message(tool_calls=[_tool_call("lookup_variant", {"gene": "NRAS", "protein_change": "G12C"})]),
        _message(content="NRAS G12C is pathogenic in ClinVar [clinvar:VCV900001] and a hotspot "
                         "[depmap:row183]. Also [clinvar:VCV12345]."),
    ])
    clients = {None: openai_client, config.GEMINI_BASE_URL: gemini_client}

    result = llm.ask("Tell me about TP53 R306*", index=index,
                     client_factory=lambda key, base_url: clients[base_url])

    assert result["provider"] == "gemini"
    assert len(openai_client.requests) == 1
    # first prompt already carries the retrieved TP53 record
    assert "depmap:row1535" in gemini_client.requests[0]["messages"][1]["content"]
    # the tool result for NRAS was fed back to the model
    assert "clinvar:VCV900001" in gemini_client.requests[1]["messages"][-1]["content"]
    assert result["tool_calls"][0]["tool"] == "lookup_variant"
    # citations are checked against what was actually retrieved
    assert result["unverified_citations"] == ["clinvar:VCV12345"]


def test_uses_openai_first(index, both_keys):
    openai_client = FakeClient(replies=[_message(content="Found [depmap:row1].")])
    result = llm.ask("PERM1 P750Q", index=index, client_factory=lambda key, base_url: openai_client)
    assert result["provider"] == "openai" and result["model"] == config.OPENAI_CHAT_MODEL
    assert result["unverified_citations"] == []


def test_gemini_only_when_no_openai_key(index, monkeypatch):
    monkeypatch.setattr(config, "OPENAI_API_KEY", None)
    monkeypatch.setattr(config, "GEMINI_API_KEY", "gm-test")
    seen = []

    def factory(key, base_url):
        seen.append(base_url)
        return FakeClient(replies=[_message(content="ok")])

    assert llm.ask("PERM1 P750Q", index=index, client_factory=factory)["provider"] == "gemini"
    assert seen == [config.GEMINI_BASE_URL]


def test_no_keys_raises(index, monkeypatch):
    monkeypatch.setattr(config, "OPENAI_API_KEY", None)
    monkeypatch.setattr(config, "GEMINI_API_KEY", None)
    with pytest.raises(llm.NoProviderError):
        llm.ask("PERM1 P750Q", index=index)
