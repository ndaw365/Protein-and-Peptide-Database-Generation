"""Build the local variant index (SQLite) from DepMap, the peptide database and ClinVar.

    python -m rag_assistant ingest                        # download ClinVar from NCBI
    python -m rag_assistant ingest --clinvar variant_summary.txt.gz
    python -m rag_assistant ingest --no-clinvar
    python -m rag_assistant ingest --gene-effect CRISPRGeneEffect.csv --depmap-model ACH-XXXXXX
    python -m rag_assistant ingest --embed                # also embed variant cards
    python -m rag_assistant embed                         # (re)start or resume embedding
"""

import csv
import gzip
import hashlib
import io
import json
import shutil
import sqlite3
import sys
import urllib.request
from pathlib import Path

from . import config
from .normalize import chrom, protein_change_key, rsid, uniprot_base
from .retrieve import FLAG_FIELDS, Index

SCHEMA = """
CREATE TABLE depmap (
    row_id INTEGER PRIMARY KEY,      -- 1-based data row in "MOLT4 mutations.csv"
    gene TEXT, change_key TEXT, protein_change TEXT, variant_info TEXT,
    uniprot TEXT, rsid TEXT, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
    data TEXT                        -- full DepMap row as JSON
);
CREATE INDEX depmap_gene_change ON depmap (gene, change_key);
CREATE INDEX depmap_rsid ON depmap (rsid);
CREATE INDEX depmap_uniprot ON depmap (uniprot);
CREATE INDEX depmap_locus ON depmap (chrom, pos, ref, alt);

CREATE TABLE peptides (
    row_id INTEGER PRIMARY KEY,      -- 1-based data row in MOLT4_peptide_data.csv
    gene TEXT, change_key TEXT, uniprot TEXT, data TEXT
);
CREATE INDEX peptides_gene_change ON peptides (gene, change_key);

CREATE TABLE dropped (
    row_id INTEGER PRIMARY KEY,      -- 1-based data row in MOLT4_dropped_variants.csv
    gene TEXT, change_key TEXT, reason TEXT
);
CREATE INDEX dropped_gene_change ON dropped (gene, change_key);

CREATE TABLE clinvar (
    variation_id TEXT, gene TEXT, change_key TEXT, rsid TEXT,
    chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
    name TEXT, significance TEXT, review_status TEXT, conditions TEXT, last_evaluated TEXT
);
CREATE INDEX clinvar_gene_change ON clinvar (gene, change_key);
CREATE INDEX clinvar_rsid ON clinvar (rsid);
CREATE INDEX clinvar_locus ON clinvar (chrom, pos, ref, alt);

CREATE TABLE detections (
    row_id INTEGER PRIMARY KEY,      -- 1-based data row in the FragPipe peptide table
    gene TEXT, change_key TEXT, peptide TEXT, start INTEGER, end_ INTEGER,
    spectral_count INTEGER, intensity REAL, match_type TEXT
);
CREATE INDEX detections_gene_change ON detections (gene, change_key);

CREATE TABLE gene_effect (gene TEXT PRIMARY KEY, model_id TEXT, effect REAL);

CREATE TABLE cards (
    card_id INTEGER PRIMARY KEY,     -- equals depmap.row_id
    gene TEXT, change_key TEXT, text TEXT, text_hash TEXT,
    -- structured fields for filter_variants; multi-valued fields are "&"-joined terms
    variant_info TEXT, vep_clin_sig TEXT, clinvar_sig TEXT, am_class TEXT, flags TEXT,
    in_peptide_db INTEGER
);

CREATE TABLE embeddings (
    card_id INTEGER PRIMARY KEY, text_hash TEXT, vector BLOB  -- float32 bytes
);

CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT);
"""


def _int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def load_depmap(conn, path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = [
            (
                i, r["Gene"], protein_change_key(r["Protein Change"]), r["Protein Change"],
                r["Variant Info"], uniprot_base(r["Uniprot ID"]), rsid(r["Dbsnp Rs ID"]),
                chrom(r["Chromosome"]), _int(r["Position"]), r["Ref Allele"], r["Alt Allele"],
                json.dumps({k: v for k, v in r.items() if v not in ("", None)}),
            )
            for i, r in enumerate(csv.DictReader(fh), start=1)
        ]
    conn.executemany("INSERT INTO depmap VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    return len(rows)


def load_peptides(conn, path):
    if not Path(path).exists():
        return 0
    with open(path, newline="", encoding="utf-8") as fh:
        rows = [
            (i, r["Gene"], protein_change_key(r["Clean_Protein_Change"]),
             uniprot_base(r.get("Uniprot.ID")), json.dumps(r))
            for i, r in enumerate(csv.DictReader(fh), start=1)
        ]
    conn.executemany("INSERT INTO peptides VALUES (?,?,?,?,?)", rows)
    return len(rows)


def load_dropped(conn, path):
    if not Path(path).exists():
        return 0
    with open(path, newline="", encoding="utf-8") as fh:
        rows = [
            (i, r["Gene"], protein_change_key(r["Protein.Change"]), r["Reason"])
            for i, r in enumerate(csv.DictReader(fh), start=1)
        ]
    conn.executemany("INSERT INTO dropped VALUES (?,?,?,?)", rows)
    return len(rows)


def download(url, dest, timeout=60):
    """Download url to dest via a .part file, so an interrupted download never leaves a
    truncated file that later runs would mistake for a complete one."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    print(f"Downloading {url} -> {dest}", file=sys.stderr)
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp, open(part, "wb") as out:
            expected = resp.headers.get("Content-Length")
            shutil.copyfileobj(resp, out, length=1 << 20)
        if expected is not None and part.stat().st_size != int(expected):
            raise OSError(f"incomplete download: got {part.stat().st_size} of {expected} bytes")
    except BaseException:
        part.unlink(missing_ok=True)
        raise
    part.replace(dest)
    return dest


def _open_clinvar(source):
    """Open a local path or URL to variant_summary.txt(.gz) as a text stream."""
    source = str(source)
    if source.startswith(("http://", "https://")):
        local = config.DATA_DIR / Path(source).name
        if not local.exists():
            download(source, local)
        source = local
    raw = gzip.open(source, "rb") if str(source).endswith(".gz") else open(source, "rb")
    return io.TextIOWrapper(raw, encoding="utf-8", newline="")


def load_clinvar(conn, source, genes):
    """Keep GRCh38 ClinVar records for genes mutated in MOLT4 (the full file is ~400 MB)."""
    try:
        return _load_clinvar_rows(conn, source, genes)
    except (EOFError, gzip.BadGzipFile, UnicodeDecodeError) as exc:
        raise SystemExit(
            f"Could not read ClinVar file {source} ({exc}). It is probably damaged: delete it "
            "and run ingest again to re-download.") from exc


def _load_clinvar_rows(conn, source, genes):
    count = 0
    with _open_clinvar(source) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        # ClinVar renamed ClinicalSignificance to GermlineClassification in 2024
        sig_col = next(
            (c for c in ("ClinicalSignificance", "GermlineClassification") if c in reader.fieldnames),
            None,
        )
        batch = []
        for r in reader:
            if r.get("Assembly") != "GRCh38" or r.get("GeneSymbol") not in genes:
                continue
            batch.append((
                r.get("VariationID"), r["GeneSymbol"], protein_change_key(r.get("Name")),
                rsid(r.get("RS# (dbSNP)")), chrom(r.get("Chromosome")), _int(r.get("PositionVCF")),
                r.get("ReferenceAlleleVCF"), r.get("AlternateAlleleVCF"), r.get("Name"),
                r.get(sig_col) if sig_col else None, r.get("ReviewStatus"),
                r.get("PhenotypeList"), r.get("LastEvaluated"),
            ))
            if len(batch) >= 5000:
                conn.executemany("INSERT INTO clinvar VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", batch)
                count += len(batch)
                batch = []
        conn.executemany("INSERT INTO clinvar VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", batch)
        count += len(batch)
    return count


def parse_variant_header(header):
    """"Fwd_spD242N|P62136-1|PPP1CA_D242N ..." -> ("PPP1CA", "D242N"); None for decoys."""
    fields = header.split()[0].split("|") if header else []
    if len(fields) < 3 or not fields[0].startswith("Fwd_sp"):
        return None
    change = fields[0].removeprefix("Fwd_sp")
    gene = fields[2].removesuffix(f"_{change}")
    return gene, protein_change_key(change)


def load_detections(conn, path):
    """FragPipe combined_peptide.tsv from a search against the variant peptide database
    (see https://github.com/ndaw365/MOLT4-Checking). Start/End are positions within the
    peptide-database entry, so they can be compared with mutation_pos_in_peptide_fwd."""
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        count_col = next(c for c in reader.fieldnames if c.endswith("Spectral Count"))
        intensity_col = next(c for c in reader.fieldnames if c.endswith("Intensity"))
        match_col = next((c for c in reader.fieldnames if c.endswith("Match Type")), None)
        rows = []
        for i, r in enumerate(reader, start=1):
            parsed = parse_variant_header(r["Protein"])
            if not parsed:
                continue
            rows.append((i, *parsed, r["Peptide Sequence"], _int(r["Start"]), _int(r["End"]),
                         _int(r[count_col]), float(r[intensity_col] or 0),
                         r[match_col] if match_col else None))
    conn.executemany("INSERT INTO detections VALUES (?,?,?,?,?,?,?,?,?)", rows)
    return len(rows)


def load_gene_effect(conn, path, model_id):
    """DepMap CRISPRGeneEffect.csv: one row per cell line, columns like "TP53 (7157)"."""
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        for row in reader:
            if row[0] == model_id:
                rows = [
                    (col.split(" (")[0], model_id, float(val))
                    for col, val in zip(header[1:], row[1:]) if val not in ("", "NA")
                ]
                conn.executemany("INSERT INTO gene_effect VALUES (?,?,?)", rows)
                return len(rows)
    raise SystemExit(f"Model {model_id} not found in {path}")


def build_cards(conn):
    """One searchable text summary per DepMap variant, used for keyword and vector search."""

    index = Index(conn)
    rows = conn.execute("SELECT row_id FROM depmap ORDER BY row_id").fetchall()
    cards = []
    for (row_id,) in rows:
        r = index.variant_record(row_id)
        cards.append((
            row_id, r["gene"], r["change_key"], r["card"], _hash(r["card"]),
            r["variant_info"], r["depmap"].get("Vep Clin Sig"),
            "&".join(c["significance"] for c in r["clinvar"] if c["significance"]) or None,
            r["depmap"].get("AM class"),
            "&".join(f for f in FLAG_FIELDS if r["depmap"].get(f) == "True") or None,
            int(r["peptide"] is not None),
        ))
    conn.executemany("INSERT INTO cards VALUES (?,?,?,?,?,?,?,?,?,?,?)", cards)
    return len(cards)


def _hash(text):
    return hashlib.sha256(text.encode()).hexdigest()


def _carry_over_embeddings(conn, old_db):
    """Keep vectors from the previous index for cards whose text has not changed."""
    try:
        conn.execute("ATTACH DATABASE ? AS old", (str(old_db),))
    except sqlite3.Error:
        return 0
    try:
        meta = dict(conn.execute("SELECT key, value FROM old.meta WHERE key LIKE 'embedding_%'"))
        if not meta.get("embedding_provider"):
            return 0
        cur = conn.execute(
            "INSERT INTO embeddings SELECT e.card_id, e.text_hash, e.vector FROM old.embeddings e"
            " JOIN cards c ON c.card_id = e.card_id AND c.text_hash = e.text_hash")
        conn.executemany("INSERT OR REPLACE INTO meta VALUES (?, ?)", meta.items())
        return cur.rowcount
    except sqlite3.Error:  # index built before embeddings were stored in it
        return 0
    finally:
        conn.commit()
        conn.execute("DETACH DATABASE old")


def embed_index(db_path=None, provider=None, batch_size=100):
    """Embed every card that has no vector yet, committing after each batch.

    Safe to interrupt or to fail on a rate limit: run it again and it resumes.
    Returns (embedded_now, total_embedded, total_cards, error_or_None).
    """
    from . import llm

    conn = sqlite3.connect(db_path or config.INDEX_DB)
    meta = dict(conn.execute("SELECT key, value FROM meta WHERE key LIKE 'embedding_%'"))
    stored = meta.get("embedding_provider")
    if provider and stored and provider != stored:
        # vectors from different models are not comparable, so start over
        conn.execute("DELETE FROM embeddings")
        stored = None
    todo = conn.execute(
        "SELECT card_id, text, text_hash FROM cards WHERE card_id NOT IN"
        " (SELECT card_id FROM embeddings) ORDER BY card_id").fetchall()
    done_now, error = 0, None
    for i in range(0, len(todo), batch_size):
        batch = todo[i:i + batch_size]
        try:
            vectors, used = llm.embed_texts([t for _, t, _ in batch], provider=provider or stored)
        except llm.NoProviderError:
            conn.close()
            raise
        except Exception as exc:  # rate limit, quota, network: keep what is done
            error = f"{type(exc).__name__}: {exc}"
            break
        if stored is None:
            conn.executemany("INSERT OR REPLACE INTO meta VALUES (?, ?)", [
                ("embedding_provider", used),
                ("embedding_model", config.GEMINI_EMBED_MODEL if used == "gemini" else config.OPENAI_EMBED_MODEL),
            ])
            stored = used
        conn.executemany("INSERT INTO embeddings VALUES (?, ?, ?)", [
            (cid, h, _to_blob(v)) for (cid, _, h), v in zip(batch, vectors)])
        conn.commit()
        done_now += len(batch)
    total = conn.execute("SELECT COUNT(*) FROM embeddings").fetchone()[0]
    cards = conn.execute("SELECT COUNT(*) FROM cards").fetchone()[0]
    conn.close()
    return done_now, total, cards, error


def _to_blob(vector):
    import numpy as np

    return np.asarray(vector, dtype="float32").tobytes()


def build_index(db_path=None, clinvar=config.CLINVAR_URL, gene_effect=None, depmap_model=None,
                embed=False, depmap_csv=None, peptide_csv=None, dropped_csv=None,
                detections=config.DETECTIONS_TSV):
    db_path = Path(db_path or config.INDEX_DB)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = db_path.with_suffix(".tmp")
    tmp.unlink(missing_ok=True)
    conn = sqlite3.connect(tmp)
    conn.executescript(SCHEMA)

    stats = {"depmap": load_depmap(conn, depmap_csv or config.DEPMAP_MUTATIONS_CSV)}
    stats["peptides"] = load_peptides(conn, peptide_csv or config.PEPTIDE_CSV)
    stats["dropped"] = load_dropped(conn, dropped_csv or config.DROPPED_CSV)
    if detections and Path(detections).exists():
        stats["detections"] = load_detections(conn, detections)
    if clinvar:
        genes = {g for (g,) in conn.execute("SELECT DISTINCT gene FROM depmap")}
        stats["clinvar"] = load_clinvar(conn, clinvar, genes)
        conn.execute("INSERT INTO meta VALUES ('clinvar_source', ?)", (str(clinvar),))
    if gene_effect:
        if not depmap_model:
            raise SystemExit("--depmap-model (MOLT4's ACH- ID from DepMap Model.csv) is required")
        stats["gene_effect"] = load_gene_effect(conn, gene_effect, depmap_model)
    stats["cards"] = build_cards(conn)
    conn.commit()
    if db_path.exists():
        stats["embeddings_reused"] = _carry_over_embeddings(conn, db_path)
    conn.close()
    # The index is saved before any embedding call, so an embedding failure cannot lose it
    tmp.replace(db_path)
    if embed:
        from .llm import NoProviderError, model_hint

        try:
            done_now, total, cards, error = embed_index(db_path)
        except NoProviderError as exc:
            stats["embedding_error"] = str(exc)
            return stats
        stats["embedded"] = f"{total}/{cards}"
        if error:
            stats["embedding_error"] = f"{error}. " + (
                model_hint(error) or "Run `python -m rag_assistant embed` to resume.")
    return stats
