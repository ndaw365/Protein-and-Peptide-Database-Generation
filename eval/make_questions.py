"""Generate eval/questions.jsonl: gold retrieval questions with their relevant records.

Three kinds, each recording "relevant", the DepMap rows a correct retrieval returns:
- single-answer: one MOLT4 variant named by 1-letter change, 3-letter HGVS, dbSNP ID or
  UniProt accession + change (also kept as "expected_source")
- gene_list: "list every variant in GENE"; relevant is all of that gene's rows
- absent: a well-known variant MOLT4 does not carry; relevant is empty and no exact match
  should be returned
"""

import csv
import json
import random
from pathlib import Path

from rag_assistant import config
from rag_assistant.normalize import AA3_TO_1, rsid

AA1_TO_3 = {v: k for k, v in AA3_TO_1.items()}
TEMPLATES = {
    "one_letter": "What is known about {gene} {change}?",
    "hgvs_3letter": "Is {gene} {hgvs3} in the peptide database?",
    "uniprot": "Show the variant {change} on UniProt {uniprot}.",
    "rsid": "Which MOLT4 variant is {rsid} and what is its clinical significance?",
    "lowercase_gene": "peptide covering {gene_lower} {change}",
}

# Well-known cancer variants that MOLT4 does not carry, in genes it has and genes it lacks
ABSENT = [("BRAF", "V600E"), ("KRAS", "G12D"), ("IDH1", "R132H"), ("JAK2", "V617F"),
          ("FLT3", "D835Y"), ("NRAS", "Q61K"), ("TP53", "R175H"), ("EGFR", "L858R")]


def to_hgvs3(change):
    ref, pos, alt = change[0], change[1:-1], change[-1]
    return f"p.{AA1_TO_3[ref]}{pos}{AA1_TO_3[alt]}"


def main(n_per_template=8, seed=7):
    with open(config.DEPMAP_MUTATIONS_CSV, newline="", encoding="utf-8") as fh:
        rows = [dict(r, row_id=i) for i, r in enumerate(csv.DictReader(fh), start=1)]
    missense = [r for r in rows if r["Variant Info"] == "missense_variant" and r["Uniprot ID"]]
    # A question is only unambiguous if its key identifies one DepMap row
    gene_change = {}
    for r in rows:
        gene_change.setdefault((r["Gene"], r["Protein Change"]), []).append(r["row_id"])
    missense = [r for r in missense if len(gene_change[(r["Gene"], r["Protein Change"])]) == 1]
    with_rsid = [r for r in rows if rsid(r["Dbsnp Rs ID"])]

    rng = random.Random(seed)
    questions = []
    for style, template in TEMPLATES.items():
        pool = with_rsid if style == "rsid" else missense
        for r in rng.sample(pool, min(n_per_template, len(pool))):
            change = r["Protein Change"].removeprefix("p.")
            questions.append({
                "style": style,
                "question": template.format(
                    gene=r["Gene"], gene_lower=r["Gene"].lower(), change=change,
                    hgvs3=to_hgvs3(change) if style == "hgvs_3letter" else "",
                    uniprot=r["Uniprot ID"].split("-")[0], rsid=rsid(r["Dbsnp Rs ID"])),
                "expected_source": f"depmap:row{r['row_id']}",
                "gene": r["Gene"],
                "protein_change": r["Protein Change"],
            })
    for q in questions:
        q["relevant"] = [q["expected_source"]]

    by_gene = {}
    for r in rows:
        by_gene.setdefault(r["Gene"], []).append(f"depmap:row{r['row_id']}")
    multi = sorted(g for g, ids in by_gene.items() if g and 2 <= len(ids) <= 9)
    for gene in rng.sample(multi, n_per_template):
        questions.append({"style": "gene_list", "question": f"List every MOLT4 variant in {gene}.",
                          "gene": gene, "relevant": by_gene[gene]})

    present = {(r["Gene"], r["Protein Change"].removeprefix("p.")) for r in rows}
    for gene, change in ABSENT:
        assert (gene, change) not in present, (gene, change)
        questions.append({"style": "absent", "question": f"Does MOLT4 carry {gene} {change}?",
                          "gene": gene, "protein_change": f"p.{change}", "relevant": []})

    out = Path(__file__).with_name("questions.jsonl")
    out.write_text("".join(json.dumps(q) + "\n" for q in questions))
    print(f"Wrote {len(questions)} questions to {out}")


if __name__ == "__main__":
    main()
