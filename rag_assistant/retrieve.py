"""Retrieval over the variant index.

Exact lookups (gene + protein change, dbSNP ID, UniProt ID) answer most variant
questions and need no model. Free-text questions fall back to BM25 keyword search
over one text card per variant, fused with vector search when embeddings exist.
"""

import json
import math
import re
import sqlite3
from collections import Counter
from pathlib import Path

from . import config
from .normalize import find_protein_changes, find_rsids, find_uniprot_ids, protein_change_key

# DepMap columns worth showing the model; the full row stays in the index
DEPMAP_FIELDS = [
    "HGNC Name", "DNA Change", "Vep Mane Select", "Allele Fraction", "Ref Count", "Alt Count",
    "Vep Impact", "Vep Clin Sig", "AM class", "AM Pathogenicity", "Revel Score", "Polyphen",
    "Sift", "Provean Prediction", "Hotspot", "Hess Driver", "Likely LOF",
    "Oncogene High Impact", "Tumor Suppressor High Impact", "CIViC Description",
    "Gnomadg AF", "Gnomade AF", "Vep Existing Variation",
]
PEPTIDE_FIELDS = [
    "peptide_sequence_fwd", "peptide_start_fwd", "peptide_end_fwd", "mutation_pos_in_peptide_fwd",
    "Peptide_Header_Fwd", "peptide_sequence_rev", "Peptide_Header_Rev",
]
_TOKEN_RE = re.compile(r"[a-z0-9_.*]+")


def embeddings_path(db_path):
    db_path = Path(db_path)
    return db_path.with_name(db_path.stem + "_embeddings.npy")


def _int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _tokens(text):
    return _TOKEN_RE.findall(text.lower())


class Index:
    def __init__(self, db=None):
        if isinstance(db, sqlite3.Connection):
            self.conn, self.db_path = db, None
        else:
            self.db_path = Path(db or config.INDEX_DB)
            if not self.db_path.exists():
                raise FileNotFoundError(
                    f"No index at {self.db_path}. Run: python -m rag_assistant ingest")
            self.conn = sqlite3.connect(self.db_path)
        self.genes = {g for (g,) in self.conn.execute("SELECT DISTINCT gene FROM depmap")}
        self._bm25 = None
        self._vectors = None

    # ---------------------------------------------------------------- records
    def variant_record(self, row_id):
        """Everything known about one DepMap variant, with a source ID for each fact."""
        (gene, change_key, protein_change, variant_info, uniprot, rs, chrom_, pos, ref, alt,
         data) = self.conn.execute(
            "SELECT gene, change_key, protein_change, variant_info, uniprot, rsid, chrom, pos,"
            " ref, alt, data FROM depmap WHERE row_id = ?", (row_id,)).fetchone()
        data = json.loads(data)
        record = {
            "source_ids": [f"depmap:row{row_id}"],
            "gene": gene,
            "protein_change": protein_change or None,
            "change_key": change_key,
            "variant_info": variant_info,
            "uniprot": uniprot,
            "rsid": rs,
            "locus_grch38": f"chr{chrom_}:{pos} {ref}>{alt}" if pos else None,
            "depmap": {k: data[k] for k in DEPMAP_FIELDS if k in data},
            "peptide": None,
            "dropped_reason": None,
            "clinvar": [],
            "ms_detection": [],
            "gene_effect": None,
        }

        if change_key:
            pep = self.conn.execute(
                "SELECT row_id, data FROM peptides WHERE gene = ? AND change_key = ?",
                (gene, change_key)).fetchone()
            if pep:
                pdata = json.loads(pep[1])
                record["peptide"] = {k: pdata[k] for k in PEPTIDE_FIELDS if pdata.get(k) not in (None, "", "NA")}
                record["source_ids"].append(f"peptide:row{pep[0]}")
            drop = self.conn.execute(
                "SELECT row_id, reason FROM dropped WHERE gene = ? AND change_key = ?",
                (gene, change_key)).fetchone()
            if drop:
                record["dropped_reason"] = drop[1]
                record["source_ids"].append(f"dropped:row{drop[0]}")

        if change_key:
            mut_pos = _int((record["peptide"] or {}).get("mutation_pos_in_peptide_fwd"))
            for det_id, peptide, start, end, spectra, intensity, match_type in self.conn.execute(
                    "SELECT row_id, peptide, start, end_, spectral_count, intensity, match_type"
                    " FROM detections WHERE gene = ? AND change_key = ? ORDER BY row_id",
                    (gene, change_key)):
                record["ms_detection"].append({
                    "peptide": peptide, "start": start, "end": end,
                    "spectral_count": spectra, "intensity": intensity, "match_type": match_type,
                    # A hit on the variant entry only supports the variant if it spans the
                    # mutated residue; otherwise it is identical to the wild-type protein.
                    "covers_variant": None if mut_pos is None else start <= mut_pos <= end,
                })
                record["source_ids"].append(f"ms:row{det_id}")

        record["clinvar"] = self._clinvar_matches(gene, change_key, rs, chrom_, pos, ref, alt)
        record["source_ids"] += [f"clinvar:VCV{c['variation_id']}" for c in record["clinvar"]]

        effect = self.conn.execute(
            "SELECT model_id, effect FROM gene_effect WHERE gene = ?", (gene,)).fetchone()
        if effect:
            record["gene_effect"] = {"model_id": effect[0], "chronos_score": round(effect[1], 3)}
            record["source_ids"].append(f"gene_effect:{gene}")

        record["card"] = self._card(record)
        return record

    def _clinvar_matches(self, gene, change_key, rs, chrom_, pos, ref, alt):
        """Most specific match wins: genomic allele, then dbSNP ID, then gene + protein change."""
        cols = "variation_id, name, significance, review_status, conditions, last_evaluated"
        queries = [
            ("same_allele", f"SELECT {cols} FROM clinvar WHERE chrom=? AND pos=? AND ref=? AND alt=?",
             (chrom_, pos, ref, alt), pos is not None),
            ("same_rsid", f"SELECT {cols} FROM clinvar WHERE rsid = ?", (rs,), rs is not None),
            ("same_protein_change", f"SELECT {cols} FROM clinvar WHERE gene = ? AND change_key = ?",
             (gene, change_key), change_key is not None),
        ]
        for match_type, sql, params, usable in queries:
            if not usable:
                continue
            rows = self.conn.execute(sql, params).fetchall()
            if rows:
                keys = ["variation_id", "name", "significance", "review_status", "conditions",
                        "last_evaluated"]
                return [dict(zip(keys, r), match=match_type) for r in rows[:5]]
        return []

    @staticmethod
    def _card(r):
        parts = [f"{r['gene']} {r['protein_change'] or ''} ({r['variant_info']})"]
        if r["depmap"].get("HGNC Name"):
            parts.append(r["depmap"]["HGNC Name"])
        for key in ("Vep Clin Sig", "AM class", "Polyphen", "Sift", "Vep Impact"):
            if r["depmap"].get(key):
                parts.append(f"{key}: {r['depmap'][key]}")
        for key in ("Hotspot", "Hess Driver", "Likely LOF", "Oncogene High Impact",
                    "Tumor Suppressor High Impact"):
            if r["depmap"].get(key) == "True":
                parts.append(key)
        if r["depmap"].get("CIViC Description"):
            parts.append(r["depmap"]["CIViC Description"])
        for c in r["clinvar"]:
            parts.append(f"ClinVar {c['significance']}: {c['conditions']}")
        if r["peptide"]:
            parts.append("in peptide database")
        if any(d["covers_variant"] for d in r["ms_detection"]):
            parts.append("variant peptide detected by mass spectrometry")
        elif r["ms_detection"]:
            parts.append("MS hits on the variant entry, none spanning the mutated residue")
        if r["dropped_reason"]:
            parts.append(f"not in peptide database: {r['dropped_reason']}")
        return ". ".join(parts)

    # ---------------------------------------------------------------- exact lookups
    def lookup(self, gene=None, protein_change=None, rsid=None, uniprot=None, limit=20):
        clauses, params = [], []
        if gene:
            clauses.append("gene = ?")
            params.append(gene.upper() if gene.upper() in self.genes else gene)
        if protein_change:
            key = protein_change_key(protein_change)
            if not key:
                return []
            clauses.append("change_key = ?")
            params.append(key)
        if rsid:
            clauses.append("rsid = ?")
            params.append(rsid.lower())
        if uniprot:
            clauses.append("uniprot = ?")
            params.append(uniprot.split("-")[0].upper())
        if not clauses:
            return []
        rows = self.conn.execute(
            f"SELECT row_id FROM depmap WHERE {' AND '.join(clauses)} ORDER BY row_id LIMIT ?",
            (*params, limit)).fetchall()
        return [self.variant_record(r) for (r,) in rows]

    def filter_variants(self, gene=None, variant_type=None, clinvar_significance=None,
                        am_class=None, flag=None, in_peptide_database=None, limit=20):
        """Structured filtering over the variant cards, e.g. all hotspot missense variants."""
        sql, params = "SELECT card_id, text FROM cards WHERE 1=1", []
        if gene:
            sql += " AND gene = ?"
            params.append(gene.upper())
        if variant_type:
            sql += " AND text LIKE ?"
            params.append(f"%({'%'.join(variant_type.split())}%")
        if clinvar_significance:
            sql += " AND (text LIKE ? OR text LIKE ?)"
            params += [f"%Vep Clin Sig: %{clinvar_significance}%", f"%ClinVar %{clinvar_significance}%"]
        if am_class:
            sql += " AND text LIKE ?"
            params.append(f"%AM class: {am_class}%")
        if flag:
            sql += " AND text LIKE ?"
            params.append(f"%. {flag}%")
        if in_peptide_database is not None:
            sql += " AND text " + ("" if in_peptide_database else "NOT ") + "LIKE '%in peptide database%'"
        rows = self.conn.execute(sql + " ORDER BY card_id LIMIT ?", (*params, limit)).fetchall()
        return [self.variant_record(r) for (r, _) in rows]

    # ---------------------------------------------------------------- search
    def _bm25_index(self):
        if self._bm25 is None:
            docs = self.conn.execute("SELECT card_id, text FROM cards ORDER BY card_id").fetchall()
            tokenised = [(cid, Counter(_tokens(text))) for cid, text in docs]
            df = Counter(t for _, tf in tokenised for t in tf)
            avg_len = sum(sum(tf.values()) for _, tf in tokenised) / max(len(tokenised), 1)
            self._bm25 = (tokenised, df, avg_len)
        return self._bm25

    def keyword_search(self, query, k=10, k1=1.2, b=0.75):
        tokenised, df, avg_len = self._bm25_index()
        n = len(tokenised)
        terms = set(_tokens(query))
        scores = []
        for cid, tf in tokenised:
            length = sum(tf.values())
            score = 0.0
            for t in terms & tf.keys():
                idf = math.log(1 + (n - df[t] + 0.5) / (df[t] + 0.5))
                score += idf * tf[t] * (k1 + 1) / (tf[t] + k1 * (1 - b + b * length / avg_len))
            if score > 0:
                scores.append((score, cid))
        return [cid for _, cid in sorted(scores, reverse=True)[:k]]

    def vector_search(self, query, k=10):
        """Cosine similarity over card embeddings; [] when no embeddings were built."""
        if self.db_path is None:
            return []
        path = embeddings_path(self.db_path)
        if not path.exists():
            return []
        import numpy as np

        from .llm import embed_texts

        if self._vectors is None:
            vectors = np.load(path)
            self._vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
            row = self.conn.execute("SELECT value FROM meta WHERE key='embedding_provider'").fetchone()
            self._embed_provider = row[0] if row else None
        try:
            (qvec,), _ = embed_texts([query], provider=self._embed_provider)
        except Exception:
            return []
        qvec = np.asarray(qvec, dtype="float32")
        sims = self._vectors @ (qvec / np.linalg.norm(qvec))
        return [int(i) + 1 for i in np.argsort(-sims)[:k]]  # card_id is 1-based

    def search(self, query, k=5):
        """Reciprocal-rank fusion of keyword and vector results."""
        fused = Counter()
        for ranking in (self.keyword_search(query, k * 4), self.vector_search(query, k * 4)):
            for rank, cid in enumerate(ranking):
                fused[cid] += 1 / (60 + rank)
        return [self.variant_record(cid) for cid, _ in fused.most_common(k)]

    # ---------------------------------------------------------------- entry point
    def extract_entities(self, question):
        words = re.findall(r"[A-Za-z0-9-]+", question)
        return {
            "genes": list(dict.fromkeys(w.upper() for w in words if w.upper() in self.genes)),
            "protein_changes": find_protein_changes(question),
            "rsids": find_rsids(question),
            "uniprot_ids": [u for u in find_uniprot_ids(question) if u not in self.genes],
        }

    def retrieve(self, question, k=5):
        """Exact lookups for any identifiers in the question, else hybrid search."""
        ents = self.extract_entities(question)
        records = []
        for rs in ents["rsids"]:
            records += self.lookup(rsid=rs)
        for change in ents["protein_changes"]:
            for gene in ents["genes"] or [None]:
                for uni in ents["uniprot_ids"] or [None]:
                    if gene or uni:
                        records += self.lookup(gene=gene, protein_change=change, uniprot=uni)
        if not records and not ents["protein_changes"]:
            for uni in ents["uniprot_ids"]:
                records += self.lookup(uniprot=uni, limit=k)
        seen, unique = set(), []
        for r in records:
            if r["source_ids"][0] not in seen:
                seen.add(r["source_ids"][0])
                unique.append(r)
        if unique:
            return {"mode": "exact", "entities": ents, "records": unique[:k]}
        return {"mode": "search", "entities": ents, "records": self.search(question, k)}
