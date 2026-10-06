
# MOLT4 Mutation Analysis and Peptide Extraction Pipeline

[![tests](https://github.com/ndaw365/Protein-and-Peptide-Database-Generation/actions/workflows/tests.yml/badge.svg)](https://github.com/ndaw365/Protein-and-Peptide-Database-Generation/actions/workflows/tests.yml)

## Overview

This pipeline processes mutation data for the MOLT4 cell line. It performs:

1. **Filtering** of missense variants with valid UniProt IDs  
2. **Fetching** canonical protein sequences via the UniProt API  
3. **Applying** amino acid mutations to sequences 
4. **Generating** forward and reverse protein sequences  
5. **Extracting** peptides around mutation sites (using K/R cleavage logic)  
6. **Saving** annotated outputs in FASTA and CSV formats  
7. **Verifying** peptide correctness through position and boundary checks

## Dependencies

Ensure the following R packages are installed:

```r
install.packages(c("tidyverse", "httr", "stringr", "readr", "knitr", "gridExtra", "RColorBrewer"))
```

## Input

**Required File:**  
`MOLT4 mutations.csv`: the MOLT4 mutation table exported from the [DepMap portal](https://depmap.org/portal/)
(cell line MOLT-4; GRCh38 coordinates with VEP, AlphaMissense, REVEL and DepMap's copy of ClinVar annotations).

**Required Columns:**
- `Uniprot ID`: UniProt accession (e.g., `Q5SV97-1`); matched case-insensitively
- `Variant Info`: Must include the string `"missense_variant"` (combined consequences such as
  `missense_variant&splice_region_variant` are kept); matched case-insensitively
- `Protein Change`: Mutation notation (e.g., `p.P750Q`); R reads it as `Protein.Change`
- `Gene`: Gene symbol (e.g., `TP53`)


## Output

| File Name                                | Description |
|------------------------------------------|-------------|
| `MOLT4_mutations_with_sequences.csv`     | Filtered entries + canonical sequences |
| `MOLT4_analysis_summary.csv`             | Variant counts at each step, from input rows to the peptide database |
| `MOLT4_mutated_protein_output.csv`       | Mutated forward/reverse protein sequences + metadata |
| `MOLT4_mutated_protein_database.fasta`   | FASTA-formatted protein sequences (mutated) |
| `MOLT4_mutated_peptide_database.fasta`   | FASTA-formatted peptides around mutations |
| `MOLT4_peptide_data.csv`                 | Peptide details and positional metadata |
| `MOLT4_dropped_variants.csv`             | Missense variants left out of the databases, with the reason |
| `MOLT4_variant_peptide_length_histogram.png` | Length distribution of the forward variant peptides |
| `uniprot_sequence_cache.csv`             | Local cache of downloaded UniProt sequences (git-ignored) |


## Usage Instructions

1. Place `MOLT4 mutations.csv` in your working directory
2. Run `Rscript MOLT4_protein_peptide_databases.R` (or pass a path:
   `Rscript MOLT4_protein_peptide_databases.R "path/to/MOLT4 mutations.csv"`), or run the
   script top to bottom in RStudio
3. Output files will be saved to your working directory. Re-runs reuse
   `uniprot_sequence_cache.csv` instead of calling UniProt again

## How the Code Works

### 1. Canonical Sequence Retrieval

- Filters for rows with valid UniProt ID and `missense_variant`
- Fetches the canonical FASTA sequence from UniProt API (retries with backoff, cached on disk)
- Adds it to the dataset

### 2. Mutation Application

- Parses mutations (e.g., `p.P750Q`)
- Replaces the amino acid at the given position in the sequence
- Validates original amino acid at the target site; variants that cannot be applied (for
  example when the Ensembl transcript and UniProt isoform disagree) are written to
  `MOLT4_dropped_variants.csv` with the reason

### 3. Reverse Sequence Generation

- Generates reverse of mutated sequence

### 4. Peptide Extraction

- Locates the second K/R residue upstream and downstream of the mutation
- Extracts peptide sequence from that window, inclusive of the K/R
- Applied to both forward and reverse sequences

### 5. Verification

- Spot-checks a reproducible random sample of 5 variants (forward and reverse peptides):
  K/R boundaries, mutation position, and that the mutated residue is in the peptide

### 6. Peptide Length Distribution

- Uses the same window as the peptide database (2nd K/R upstream to 2nd K/R downstream),
  skipping the 104 variants with fewer than 2 K/R on one side, which leaves 1765 peptides
- Counts forward peptides only: each reverse (decoy) peptide has the same length as its
  forward peptide
- Bins: below 7 aa (11, 0.6%), 7–50 aa (1349, 76.4%), above 50 aa (405, 22.9%). The plot's
  x-axis stops at 100 aa, so the longest peptides are counted but not drawn

## FASTA Header Example (forward strand) for protein database:
>Fwd_spP750Q|Q5SV97-1|PERM1_P750Q OS=Homo sapiens GN=PERM1 (Sequence)

## FASTA Header Example (forward strand) for peptide database:
>Fwd_spP750Q|Q5SV97-1|PERM1_P750Q OS=Homo sapiens GN=PERM1 (Truncated Sequence)

## Downstream proteomics check

The peptide and protein databases were searched against MOLT4 bottom-up proteomics data
with FragPipe, and the results were compared in
[ndaw365/MOLT4-Checking](https://github.com/ndaw365/MOLT4-Checking). The peptide hits from that
search (`matched_peptides_to_tryptic_peptides.tsv`) are copied here as
`data_sources/MOLT4_detected_variant_peptides.tsv` so the assistant below can report them.

That search used the database as it was before the filter fix (1810 variants), so the 59
variants recovered by the fix have not been searched yet.

Of the 41 target peptide hits, only 7 span the mutated residue (EVL D20N, COPS7B A224T,
TUBA4A E77D, CDC45 E259K, DUSP7 R102C, BRD4 P1131A, P2RX4 Y378F). The other 34 match the
variant entry but are identical to the wild-type protein, so they are not evidence for the
variant. The assistant reports this as `covers_variant`.

## Variant lookup assistant (RAG)

`rag_assistant/` is a retrieval-augmented assistant for questions about MOLT4 variants. It
indexes the following sources in a local SQLite database, and every fact it returns carries
a source ID:

| Source | Source IDs |
|--------|------------|
| DepMap MOLT4 mutations (`MOLT4 mutations.csv`) | `depmap:rowN` |
| Variant peptide database (`MOLT4_peptide_data.csv`) | `peptide:rowN` |
| Dropped variants (`MOLT4_dropped_variants.csv`) | `dropped:rowN` |
| ClinVar `variant_summary.txt.gz` (GRCh38, MOLT4 genes only) | `clinvar:VCV<VariationID>` |
| FragPipe hits on the variant peptide database | `ms:rowN` |
| DepMap CRISPR gene effect (optional) | `gene_effect:GENE` |

How retrieval works:

- **Exact lookup.** Gene + protein change in 1- or 3-letter HGVS (`P750Q`, `p.Pro750Gln`,
  `R306Ter`, `K267fs`), dbSNP rsID, or UniProt accession. This needs no model. A question
  that names a gene but no change ("all variants in KMT2D") lists that gene's variants.
- **Filters** (ClinVar significance, AlphaMissense class, variant type, DepMap flags such as
  Hotspot) match whole terms: `pathogenic` matches `pathogenic&likely_pathogenic` but not
  `likely_pathogenic`, and never matches a different annotation field.
- **ClinVar matching.** A variant is linked to ClinVar by the same GRCh38 allele first,
  then the same rsID, then the same gene + protein change. The match type is reported.
- **Free-text questions** use BM25 keyword search over one summary per variant, combined
  with vector search when embeddings have been built. Embeddings are stored in the index
  and saved in batches, so a rate limit or network error never loses the index or the
  embeddings already made: run `python -m rag_assistant embed` again to continue.
  Rebuilding the index keeps the embeddings of variants whose summary did not change.
- **Answers.** The LLM gets the retrieved records plus lookup tools. It must cite source IDs,
  and any citation that was not actually retrieved is flagged. OpenAI is used first; if
  there is no `OPENAI_API_KEY` or the call fails, the same request goes to Gemini through
  its OpenAI-compatible API.

```bash
pip install -r requirements.txt
cp .env.example .env   # then put your OPENAI_API_KEY and/or GEMINI_API_KEY in .env
                       # (.env is git-ignored and loaded automatically; keys exported
                       #  in the shell take priority)

python -m rag_assistant ingest                       # downloads ClinVar from NCBI (~400 MB)
python -m rag_assistant ingest --clinvar variant_summary.txt.gz --embed
python -m rag_assistant embed                        # resume embedding after a rate limit
# optional: --gene-effect CRISPRGeneEffect.csv --depmap-model <MOLT4 ACH- ID from DepMap Model.csv>

python -m rag_assistant lookup NRAS p.Gly12Cys       # exact lookup, no LLM
python -m rag_assistant lookup --rsid rs121913250
python -m rag_assistant search "pathogenic tumor suppressor stop gained"
python -m rag_assistant ask "Was the EVL D20N variant peptide detected by mass spec?"
python -m rag_assistant ask "..." --no-llm           # show retrieved records only
python -m rag_assistant ask "..." --provider gemini  # force one provider
```

### Tests and evaluation

```bash
python -m pytest                       # synthetic ClinVar fixture, mocked LLM clients, no keys
                                       # (also run by GitHub Actions on every push)
python -m eval.make_questions          # regenerate eval/questions.jsonl (40 questions, 5 phrasings)
python eval/run.py                     # retrieval accuracy and latency vs. scanning the CSV
python eval/run.py --llm               # also check that LLM answers cite the right record
```

On the committed index inputs, `eval/run.py` gives 100% top-1 retrieval on the 40 questions
(median about 0.3 ms per query). A plain scan of the CSV gets 20%: it only finds the exact
1-letter phrasing, not 3-letter HGVS, rsID, UniProt or lower-case gene questions, and takes
about 35 ms per query. The questions are generated from the same data, so treat this as a
check of the lookup logic, not a benchmark of open-ended questions.

