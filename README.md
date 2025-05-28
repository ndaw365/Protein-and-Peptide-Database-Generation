
# MOLT4 Mutation Analysis and Peptide Extraction Pipeline

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
install.packages(c("tidyverse", "httr", "stringr", "readr", "knitr"))
```

## Input

**Required File:**  
`MOLT4 mutations.csv`

**Required Columns (case-insensitive):**
- `Uniprot ID`: UniProt accession (e.g., P12345)
- `Variant Info`: Must include the string `"missense_variant"`
- `Protein.Change`: Mutation notation (e.g., `p.P750Q`)
- `Gene`: Gene symbol (e.g., `TP53`)


## Output

| File Name                                | Description |
|------------------------------------------|-------------|
| `MOLT4_mutations_with_sequences.csv`     | Filtered entries + canonical sequences |
| `MOLT4_analysis_summary.csv`             | Summary of processing/filtering stats |
| `MOLT4_mutated_protein_output.csv`       | Mutated forward/reverse protein sequences + metadata |
| `MOLT4_mutated_protein_database.fasta`   | FASTA-formatted protein sequences (mutated) |
| `MOLT4_mutated_peptide_database.fasta`   | FASTA-formatted peptides around mutations |
| `MOLT4_peptide_data.csv`                 | Peptide details and positional metadata |


## Usage Instructions

1. Place `MOLT4 mutations.csv` in your working directory
2. Open the R script or R Markdown file
3. Run all code blocks from top to bottom (recommended: RStudio)
4. Output files will be saved to your working directory

## How the Code Works

### 1. Canonical Sequence Retrieval

- Filters for rows with valid UniProt ID and `missense_variant`
- Fetches the canonical FASTA sequence from UniProt API
- Adds it to the dataset

### 2. Mutation Application

- Parses mutations (e.g., `p.P750Q`)
- Replaces the amino acid at the given position in the sequence
- Validates original amino acid at the target site

### 3. Reverse Sequence Generation

- Generates reverse of mutated sequence

### 4. Peptide Extraction

- Locates the second K/R residue upstream and downstream of the mutation
- Extracts peptide sequence from that window, inclusive of the K/R
- Applied to both forward and reverse sequences

### 5. Verification

- Checks boundary correctness of extracted peptides
- Validates mutation position
- Ensures the mutation is reflected in the peptide

## FASTA Header Example (forward strand) for protein database:
>Fwd_spP750Q|Q5SV97-1|PERM1_P750Q OS=Homo sapiens GN=PERM1P750Q (Sequence)

## FASTA Header Example (forward strand) for peptide database:
>Fwd_spP750Q|Q5SV97-1|PERM1_P750Q OS=Homo sapiens GN=PERM1 (Truncated Sequence)












