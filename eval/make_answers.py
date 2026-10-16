"""Generate answer-eval sets from the index, so their expected values come from the data.

    python -m eval.make_answers            # after python -m rag_assistant ingest

Writes:
- eval/answers_generated.jsonl: fact cases (peptide sequence, peptide start/end,
  AlphaMissense class, UniProt accession) for randomly chosen MOLT4 missense variants.
  Run them with: python eval/run.py --llm --answers eval/answers_generated.jsonl
- eval/judge_calibration.jsonl: evidence + answer pairs labelled faithful or not, with
  planted errors, to measure how well the faithfulness judge catches them
  (python eval/run.py --calibrate-judge).
"""

import argparse
import json
import random
from pathlib import Path

from rag_assistant import config
from rag_assistant.llm import _compact
from rag_assistant.retrieve import Index

EVAL_DIR = Path(__file__).resolve().parent

TEMPLATES = {
    "peptide": "What peptide covers {gene} {change} in the peptide database?",
    "positions": "Where does the peptide for {gene} {change} start and end in the protein?",
    "am_class": "What is the AlphaMissense class of {gene} {change}?",
    "uniprot": "Which UniProt accession is {gene} {change} annotated on?",
}


def _case(kind, record):
    gene, change = record["gene"], record["protein_change"].removeprefix("p.")
    pep = record["peptide"]
    depmap_id = record["source_ids"][0]
    peptide_id = next((s for s in record["source_ids"] if s.startswith("peptide:")), None)
    case = {"id": f"gen-{kind}-{gene}-{change}", "category": "fact",
            "question": TEMPLATES[kind].format(gene=gene, change=change), "relevant": [depmap_id]}
    if kind == "peptide":
        case.update(must_include=[[pep["peptide_sequence_fwd"]]], expected_source=peptide_id)
    elif kind == "positions":
        case.update(must_include=[[pep["peptide_start_fwd"]], [pep["peptide_end_fwd"]]],
                    expected_source=peptide_id)
    elif kind == "am_class":
        am = record["depmap"]["AM class"]
        case.update(must_include=[[am, am.replace("_", " ")]], expected_source=depmap_id)
    else:
        case.update(must_include=[[record["uniprot"]]], expected_source=depmap_id)
    return case


def generated_cases(index, n_per_template=8, seed=11):
    rng = random.Random(seed)
    rows = index.conn.execute(
        "SELECT row_id, gene, protein_change FROM depmap WHERE variant_info = 'missense_variant'"
        " ORDER BY row_id").fetchall()
    records = []
    for row_id, gene, change in rows:
        matches = index.lookup(gene=gene, protein_change=change)
        if len(matches) == 1 and matches[0]["peptide"]:  # unambiguous and in the database
            records.append(matches[0])
    cases = []
    for kind in TEMPLATES:
        pool = [r for r in records if kind != "am_class" or r["depmap"].get("AM class")]
        cases += [_case(kind, r) for r in rng.sample(pool, n_per_template)]
    return cases


# (gene, protein change or None for all of the gene's variants, answer, faithful, planted error)
CALIBRATION = [
    ("PERM1", "P750Q", "The peptide covering PERM1 P750Q is "
     "RGPVPSFAFSQNDMCLVFVAFATWAVRTSDQHTPDAWKTALLANVGTISAIR [peptide:row1].", True, None),
    ("PERM1", "P750Q", "The peptide covering PERM1 P750Q is "
     "RGPVPSFAFSQNDMCLVFVAFATWAVRTSDQHTPDAWKTALLANVGTISAIK [peptide:row1].", False,
     "wrong peptide residue"),
    ("EVL", "D20N", "EVL D20N was detected by mass spectrometry: the peptide ASVMVYDNTSK spans "
     "the mutated residue [ms:row3].", True, None),
    ("PPP1CA", "D242N", "PPP1CA D242N was detected by mass spectrometry, confirming the variant "
     "protein is expressed [ms:row1].", False, "MS hit misread as covering the variant"),
    ("PPP1CA", "D242N", "PPP1CA D242N has an MS hit (AHQVVEDGYEFFAK), but that peptide does not "
     "span the mutated residue, so it is not evidence for the variant [ms:row1].", True, None),
    ("NRAS", "G12C", "NRAS G12C is a hotspot, and DepMap's ClinVar annotation is "
     "pathogenic/likely pathogenic [depmap:row183].", True, None),
    ("NRAS", "G12C", "NRAS G12C is a hotspot [depmap:row183] and can be treated with "
     "sotorasib.", False, "outside knowledge (drug)"),
    ("HELZ2", "R563L", "HELZ2 R563L is not in the peptide database because UniProt has A, not R, "
     "at position 563 [dropped:row1].", True, None),
    ("HELZ2", "R563L", "HELZ2 R563L is not in the peptide database because it is a frameshift "
     "[dropped:row1].", False, "wrong drop reason"),
    ("KMT2D", None, "MOLT4 has 9 variants in KMT2D [depmap:row795].", True, None),
    ("KMT2D", None, "MOLT4 has 12 variants in KMT2D [depmap:row795].", False, "wrong count"),
    ("STK11", "M125I", "STK11 M125I is in the peptide database [peptide:row953]; its AlphaMissense "
     "class is likely pathogenic and DepMap's ClinVar annotation is uncertain significance "
     "[depmap:row1916].", True, None),
    ("STK11", "M125I", "STK11 M125I is in the peptide database [peptide:row953] and ClinVar "
     "classifies it as pathogenic [depmap:row1916].", False, "wrong ClinVar call"),
]


def calibration_items(index):
    items = []
    for n, (gene, change, answer, faithful, error) in enumerate(CALIBRATION, start=1):
        records = index.lookup(gene=gene, protein_change=change, limit=50)
        evidence = [_compact(r) for r in records]
        for r in evidence:  # keep items independent of the ClinVar release
            r.pop("clinvar", None)
            r["source_ids"] = [s for s in r["source_ids"] if not s.startswith("clinvar:")]
        question = f"Tell me about {gene} {change}" if change else f"List the {gene} variants"
        items.append({"id": f"cal-{n}", "question": question, "answer": answer,
                      "evidence": evidence, "faithful": faithful, "planted_error": error})
    return items


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default=str(config.INDEX_DB))
    args = parser.parse_args()
    index = Index(args.db)
    for name, rows in (("answers_generated.jsonl", generated_cases(index)),
                       ("judge_calibration.jsonl", calibration_items(index))):
        (EVAL_DIR / name).write_text("".join(json.dumps(r) + "\n" for r in rows))
        print(f"Wrote {len(rows)} items to eval/{name}")


if __name__ == "__main__":
    main()
