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
