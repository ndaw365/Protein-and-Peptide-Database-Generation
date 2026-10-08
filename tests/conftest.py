from pathlib import Path

import pytest

from rag_assistant.ingest import build_index
from rag_assistant.retrieve import Index

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="session")
def index(tmp_path_factory):
    """Index built from the real MOLT4 files plus the synthetic ClinVar fixture."""
    db = tmp_path_factory.mktemp("index") / "variants.sqlite"
    build_index(db_path=db, clinvar=FIXTURES / "clinvar_sample.tsv")
    return Index(db)


@pytest.fixture(autouse=True)
def no_real_llm_calls(monkeypatch):
    """Tests must never reach a real provider, even with keys in .env."""
    from rag_assistant import llm

    def blocked(*_, **__):
        raise AssertionError("a test tried to create a real LLM client")

    monkeypatch.setattr(llm, "_client", blocked)
