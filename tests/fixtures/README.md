`clinvar_sample.tsv` is **synthetic test data** in the column layout of ClinVar's
`variant_summary.txt` (2024+ format, with `GermlineClassification`). VariationIDs
900001–900004 and the "SYNTHETIC test condition" phenotypes are made up: they exist only to
exercise the three ClinVar match paths (same allele, same rsID, same protein change) and the
GRCh38 / MOLT4-gene filters. Do not read them as real ClinVar classifications.
