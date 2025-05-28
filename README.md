
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

- **`MOLT4 mutations.csv`**  
  Required columns (case-insensitive):
  - `Uniprot ID`: valid UniProt accession (e.g., P12345)
  - `Variant Info`: must contain `"missense_variant"`
  - `Protein.Change`: mutation in format like `p.P750Q`
  - `Gene`: gene name (e.g., TP53)

## Output

| File Name                                | Description |
|------------------------------------------|-------------|
| `MOLT4_mutations_with_sequences.csv`     | Filtered entries + canonical sequences |
| `MOLT4_analysis_summary.csv`             | Summary of processing/filtering stats |
| `MOLT4_mutated_protein_output.csv`       | Mutated forward/reverse protein sequences + metadata |
| `MOLT4_mutated_protein_database.fasta`   | FASTA-formatted protein sequences (mutated) |
| `MOLT4_mutated_peptide_database_CORRECTED.fasta` | FASTA-formatted peptides around mutations |
| `MOLT4_peptide_data.csv`                 | Peptide details and positional metadata |













