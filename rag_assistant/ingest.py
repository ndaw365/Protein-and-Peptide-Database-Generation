"""Build the local variant index (SQLite) from DepMap, the peptide database and ClinVar.

    python -m rag_assistant ingest                        # download ClinVar from NCBI
    python -m rag_assistant ingest --clinvar variant_summary.txt.gz
    python -m rag_assistant ingest --no-clinvar
    python -m rag_assistant ingest --gene-effect CRISPRGeneEffect.csv --depmap-model ACH-XXXXXX
    python -m rag_assistant ingest --embed                # also embed variant cards
"""

import csv
import gzip
import io
import json
import shutil
import sqlite3
import sys
import urllib.request
from pathlib import Path

from . import config
from .normalize import chrom, protein_change_key, rsid, uniprot_base
from .retrieve import Index, embeddings_path

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
    gene TEXT, change_key TEXT, text TEXT
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


def _open_clinvar(source):
    """Open a local path or URL to variant_summary.txt(.gz) as a text stream."""
    source = str(source)
    if source.startswith(("http://", "https://")):
        config.DATA_DIR.mkdir(parents=True, exist_ok=True)
        local = config.DATA_DIR / Path(source).name
        if not local.exists():
            print(f"Downloading {source} -> {local}", file=sys.stderr)
            with urllib.request.urlopen(source) as resp, open(local, "wb") as out:
                shutil.copyfileobj(resp, out)
        source = local
    raw = gzip.open(source, "rb") if str(source).endswith(".gz") else open(source, "rb")
    return io.TextIOWrapper(raw, encoding="utf-8", newline="")


def load_clinvar(conn, source, genes):
    """Keep GRCh38 ClinVar records for genes mutated in MOLT4 (the full file is ~400 MB)."""
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
        record = index.variant_record(row_id)
        cards.append((row_id, record["gene"], record["change_key"], record["card"]))
    conn.executemany("INSERT INTO cards VALUES (?,?,?,?)", cards)
    return len(cards)


def embed_cards(conn, out_path):
    import numpy as np

    from .llm import embed_texts

    texts = [t for (t,) in conn.execute("SELECT text FROM cards ORDER BY card_id")]
    vectors, provider = embed_texts(texts)
    np.save(out_path, np.asarray(vectors, dtype="float32"))
    conn.execute("INSERT OR REPLACE INTO meta VALUES ('embedding_provider', ?)", (provider,))
    return len(texts), provider


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
    vectors = embeddings_path(db_path)
    vectors.unlink(missing_ok=True)  # stale vectors would no longer line up with the cards
    if embed:
        stats["embedded"], stats["embedding_provider"] = embed_cards(conn, vectors)
        conn.commit()
    conn.close()
    tmp.replace(db_path)
    return stats
