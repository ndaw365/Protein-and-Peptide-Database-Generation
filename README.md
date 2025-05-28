# MOLT4 Mutation Analysis and Peptide Extraction Pipeline

## Overview

This pipeline processes a CSV file containing mutation data for the MOLT4 cell line. It performs the following main tasks:

1. Filters for missense variants with valid UniProt IDs

2. Fetches canonical protein sequences from UniProt API

3. Applies amino acid mutations to the sequences

4. Generates forward and reverse sequences

5. Extracts peptides flanked by two cleavage (K/R) sites around mutations

6. Outputs annotated FASTA and CSV files

7. Performs validation of peptide extraction

## Inputs

**MOLT4 mutations.csv:** Input CSV file with mutation data. The following columns must be present:

Uniprot ID: Identifies the protein sequence in UniProt (case-insensitive column name match).

Variant Info: Must include the string missense_variant to be considered for processing.

Protein.Change: Denotes the amino acid change (e.g., p.P750Q).

Gene: The gene symbol related to the protein mutation.

## Outputs

The following files will be generated:

**MOLT4_mutations_with_sequences.csv:** Contains the original input with appended canonical sequences.

**MOLT4_analysis_summary.csv:** Summary statistics about how many entries were filtered, matched, and retrieved.

**MOLT4_mutated_protein_output.csv:** Sequences with mutations applied, forward and reverse sequences, and sequence headers.

**MOLT4_mutated_protein_database.fasta:** FASTA file with both forward and reverse mutated protein sequences.

**MOLT4_mutated_peptide_database_CORRECTED.fasta:** Peptides extracted around mutation sites, flanked by the second K/R residue on each side.

**MOLT4_peptide_data.csv:** Summary of peptide extraction, positions, headers, and mutation validation.



# 🧬 MOLT4 Mutation Analysis and Peptide Extraction Pipeline

## 📌 Overview

This pipeline processes mutation data for the MOLT4 cell line. It performs:

1. **Filtering** of missense variants with valid UniProt IDs  
2. **Fetching** canonical protein sequences via the UniProt API  
3. **Applying** amino acid mutations to sequences  
4. **Generating** forward and reverse protein sequences  
5. **Extracting** peptides around mutation sites (using K/R cleavage logic)  
6. **Saving** annotated outputs in FASTA and CSV formats  
7. **Verifying** peptide correctness through position and boundary checks

---

## Input

- **`MOLT4 mutations.csv`**  
  Required columns (case-insensitive):
  - `Uniprot ID`: valid UniProt accession (e.g., P12345)
  - `Variant Info`: must contain `"missense_variant"`
  - `Protein.Change`: mutation in format like `p.P750Q`
  - `Gene`: gene name (e.g., TP53)

---

## Output

| File Name                                | Description |
|------------------------------------------|-------------|
| `MOLT4_mutations_with_sequences.csv`     | Filtered entries + canonical sequences |
| `MOLT4_analysis_summary.csv`             | Summary of processing/filtering stats |
| `MOLT4_mutated_protein_output.csv`       | Mutated forward/reverse protein sequences + metadata |
| `MOLT4_mutated_protein_database.fasta`   | FASTA-formatted protein sequences (mutated) |
| `MOLT4_mutated_peptide_database_CORRECTED.fasta` | FASTA-formatted peptides around mutations |
| `MOLT4_peptide_data.csv`                 | Peptide details and positional metadata |

---

## Dependencies

Ensure the following R packages are installed:

```r
install.packages(c("tidyverse", "httr", "stringr", "readr", "knitr"))












