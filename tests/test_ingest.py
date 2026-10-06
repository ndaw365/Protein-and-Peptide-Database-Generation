"""Index building: ClinVar download, damaged files, and resumable embeddings."""

import gzip
import io
import sqlite3

import numpy as np
import pytest

from rag_assistant import ingest, llm
from rag_assistant.retrieve import Index

from .conftest import FIXTURES


class FakeResponse(io.BytesIO):
    def __init__(self, body, content_length=None, fail_after=None):
        super().__init__(body)
        self.headers = {"Content-Length": str(content_length if content_length is not None else len(body))}
        self.fail_after = fail_after

    def read(self, size=-1):
        if self.fail_after is not None and self.tell() >= self.fail_after:
            raise ConnectionResetError("connection dropped")
        return super().read(size if self.fail_after is None else min(size, self.fail_after))

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_download_writes_file_only_when_complete(tmp_path, monkeypatch):
    body = b"x" * 5000
    monkeypatch.setattr(ingest.urllib.request, "urlopen", lambda url, timeout: FakeResponse(body))
    dest = ingest.download("https://example.org/f.gz", tmp_path / "f.gz")
    assert dest.read_bytes() == body
    assert not (tmp_path / "f.gz.part").exists()


@pytest.mark.parametrize("response", [
    lambda: FakeResponse(b"x" * 100, content_length=5000),     # server closed early
    lambda: FakeResponse(b"x" * 5000, fail_after=1000),        # connection dropped mid-way
])
def test_interrupted_download_leaves_nothing_behind(tmp_path, monkeypatch, response):
    monkeypatch.setattr(ingest.urllib.request, "urlopen", lambda url, timeout: response())
    with pytest.raises(OSError):
        ingest.download("https://example.org/f.gz", tmp_path / "f.gz")
    assert list(tmp_path.iterdir()) == []


def test_damaged_clinvar_file_gives_a_clear_error(tmp_path):
    good = gzip.compress((FIXTURES / "clinvar_sample.tsv").read_bytes())
    damaged = tmp_path / "variant_summary.txt.gz"
    damaged.write_bytes(good[: len(good) // 2])
    with pytest.raises(SystemExit, match="delete it"):
        ingest.build_index(db_path=tmp_path / "v.sqlite", clinvar=damaged)


def _fake_vector(text):
    """Deterministic stand-in for an embedding: letter counts."""
    vec = np.zeros(26, dtype="float32")
    for ch in text.lower():
        if "a" <= ch <= "z":
            vec[ord(ch) - 97] += 1
    return vec + 1e-3


def fake_embedder(fail_after_calls=None):
    calls = {"n": 0}

    def embed_texts(texts, provider=None, **_):
        calls["n"] += 1
        if fail_after_calls is not None and calls["n"] > fail_after_calls:
            raise RuntimeError("All LLM providers failed: gemini: RateLimitError")
        return [_fake_vector(t) for t in texts], provider or "gemini"

    return embed_texts


def test_embedding_failure_keeps_the_index_and_resumes(tmp_path, monkeypatch):
    db = tmp_path / "v.sqlite"
    monkeypatch.setattr(llm, "embed_texts", fake_embedder(fail_after_calls=3))
    stats = ingest.build_index(db_path=db, clinvar=None, embed=True)

    # the index was saved even though embedding stopped part-way
    assert stats["embedded"] == "300/3826"
    assert "resume" in stats["embedding_error"]
    assert Index(db).lookup(gene="NRAS", protein_change="G12C")

    monkeypatch.setattr(llm, "embed_texts", fake_embedder())
    done_now, total, cards, error = ingest.embed_index(db)
    assert (done_now, total, cards, error) == (3526, 3826, 3826, None)


def test_rebuild_reuses_embeddings_for_unchanged_cards(tmp_path, monkeypatch):
    db = tmp_path / "v.sqlite"
    monkeypatch.setattr(llm, "embed_texts", fake_embedder())
    ingest.build_index(db_path=db, clinvar=None, embed=True)

    # adding ClinVar changes only the cards of the 3 variants it matches
    monkeypatch.setattr(llm, "embed_texts", fake_embedder(fail_after_calls=0))
    stats = ingest.build_index(db_path=db, clinvar=FIXTURES / "clinvar_sample.tsv")
    assert stats["embeddings_reused"] == 3826 - 3
    conn = sqlite3.connect(db)
    assert dict(conn.execute("SELECT key, value FROM meta WHERE key = 'embedding_provider'")) == {
        "embedding_provider": "gemini"}


def test_vector_search_uses_stored_embeddings(tmp_path, monkeypatch):
    db = tmp_path / "v.sqlite"
    monkeypatch.setattr(llm, "embed_texts", fake_embedder())
    ingest.build_index(db_path=db, clinvar=None, embed=True)
    index = Index(db)
    target = index.conn.execute("SELECT card_id, text FROM cards WHERE card_id = 183").fetchone()
    assert index.vector_search(target[1], k=1) == [183]


def test_vector_search_without_embeddings_returns_nothing(index):
    assert index.vector_search("anything") == []
