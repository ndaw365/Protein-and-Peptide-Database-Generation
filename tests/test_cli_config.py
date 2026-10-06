import os

import pytest

from rag_assistant import cli, config

from .conftest import FIXTURES


@pytest.fixture(scope="module")
def db(tmp_path_factory):
    path = tmp_path_factory.mktemp("cli") / "v.sqlite"
    assert cli.main(["--db", str(path), "ingest", "--clinvar", str(FIXTURES / "clinvar_sample.tsv")]) == 0
    return str(path)


def test_lookup_command(db, capsys):
    assert cli.main(["--db", db, "lookup", "NRAS", "p.Gly12Cys"]) == 0
    out = capsys.readouterr().out
    assert "depmap:row183" in out and "clinvar:VCV900001" in out


def test_lookup_by_rsid(db, capsys):
    cli.main(["--db", db, "lookup", "--rsid", "rs121913344"])
    assert "TP53 p.R306Ter" in capsys.readouterr().out


def test_search_command(db, capsys):
    assert cli.main(["--db", db, "search", "tumor suppressor stop gained", "-k", "3"]) == 0
    assert capsys.readouterr().out.count("\n== ") == 3


def test_ask_no_llm_reports_partial_gene_lists(db, capsys):
    assert cli.main(["--db", db, "ask", "List every KMT2D variant", "--no-llm"]) == 0
    out = capsys.readouterr().out
    assert "Retrieval mode: exact" in out and "KMT2D" in out


def test_ask_without_keys_exits_cleanly(db, capsys, monkeypatch):
    monkeypatch.setattr(config, "OPENAI_API_KEY", None)
    monkeypatch.setattr(config, "GEMINI_API_KEY", None)
    assert cli.main(["--db", db, "ask", "NRAS G12C"]) == 2
    assert "GEMINI_API_KEY" in capsys.readouterr().err


def test_missing_index_message(tmp_path):
    with pytest.raises(FileNotFoundError, match="rag_assistant ingest"):
        cli.main(["--db", str(tmp_path / "none.sqlite"), "lookup", "NRAS", "G12C"])


def test_load_dotenv(tmp_path, monkeypatch):
    env = tmp_path / ".env"
    env.write_text(
        "\ufeff# comment, after a byte-order mark as Windows Notepad writes\n"
        "GEMINI_API_KEY=\"gm-from-file\"\n"
        "export RAG_TEST_EXPORTED='quoted value'\n"
        "OPENAI_API_KEY=\n"
        "RAG_TEST_ALREADY_SET=from-file\n"
        "not a variable line\n")
    for key in ("GEMINI_API_KEY", "RAG_TEST_EXPORTED", "OPENAI_API_KEY"):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setenv("RAG_TEST_ALREADY_SET", "from-shell")

    config.load_dotenv(env)

    assert os.environ["GEMINI_API_KEY"] == "gm-from-file"
    assert os.environ["RAG_TEST_EXPORTED"] == "quoted value"
    assert "OPENAI_API_KEY" not in os.environ          # empty values are skipped
    assert os.environ["RAG_TEST_ALREADY_SET"] == "from-shell"  # the shell wins
    config.load_dotenv(tmp_path / "missing.env")        # no file: no error


def test_embed_without_keys_says_so(db, capsys, monkeypatch):
    monkeypatch.setattr(config, "OPENAI_API_KEY", None)
    monkeypatch.setattr(config, "GEMINI_API_KEY", None)
    assert cli.main(["--db", db, "embed"]) == 2
    err = capsys.readouterr().err
    assert "No API key" in err and "resume" not in err
