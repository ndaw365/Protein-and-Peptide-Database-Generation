"""Paths and model settings. Every value can be overridden with an environment variable."""

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def load_dotenv(path):
    """Read KEY=value lines from a .env file into os.environ.

    Variables already set in the environment win, so a key exported in the shell
    overrides the file. Blank lines and # comments are ignored; values may be quoted.
    """
    path = Path(path)
    if not path.is_file():
        return
    # utf-8-sig drops the byte-order mark Windows Notepad may add
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.removeprefix("export ").split("=", 1)
        key, value = key.strip(), value.strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in "'\"":
            value = value[1:-1]
        if key and value and key not in os.environ:
            os.environ[key] = value


load_dotenv(REPO_ROOT / ".env")

# Pipeline inputs and outputs (written by MOLT4_protein_peptide_databases.R)
DEPMAP_MUTATIONS_CSV = REPO_ROOT / "MOLT4 mutations.csv"
PEPTIDE_CSV = REPO_ROOT / "MOLT4_peptide_data.csv"
DROPPED_CSV = REPO_ROOT / "MOLT4_dropped_variants.csv"
# FragPipe peptide hits against the variant peptide database (from ndaw365/MOLT4-Checking)
DETECTIONS_TSV = REPO_ROOT / "data_sources" / "MOLT4_detected_variant_peptides.tsv"

# Raw downloads and built indexes (git-ignored)
DATA_DIR = Path(os.environ.get("RAG_DATA_DIR", REPO_ROOT / "data"))
INDEX_DB = DATA_DIR / "variants.sqlite"

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

# LLM providers: OpenAI first, Gemini (via its OpenAI-compatible endpoint) as fallback
OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY")
OPENAI_CHAT_MODEL = os.environ.get("OPENAI_CHAT_MODEL", "gpt-4o-mini")
OPENAI_EMBED_MODEL = os.environ.get("OPENAI_EMBED_MODEL", "text-embedding-3-small")

GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY") or os.environ.get("GOOGLE_API_KEY")
GEMINI_BASE_URL = os.environ.get(
    "GEMINI_BASE_URL", "https://generativelanguage.googleapis.com/v1beta/openai/"
)
GEMINI_CHAT_MODEL = os.environ.get("GEMINI_CHAT_MODEL", "gemini-2.5-flash")
GEMINI_EMBED_MODEL = os.environ.get("GEMINI_EMBED_MODEL", "gemini-embedding-001")
