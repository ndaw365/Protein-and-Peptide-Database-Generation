# MOLT4 Variant Protein & Peptide Databases

[![tests](https://github.com/ndaw365/Protein-and-Peptide-Database-Generation/actions/workflows/tests.yml/badge.svg)](https://github.com/ndaw365/Protein-and-Peptide-Database-Generation/actions/workflows/tests.yml)

Builds protein and peptide search databases from MOLT4 missense variants. Also includes an
assistant that answers questions about those variants, citing DepMap, ClinVar and mass spec
results.

## Workflow

```
Stage 1  MOLT4 mutations.csv (DepMap) ── R pipeline ──► protein + peptide FASTA (variant + decoy)
                                                          └─ FragPipe search (MOLT4-Checking) ──► detected peptides
Stage 2  DepMap + peptide DB + ClinVar + detected peptides ── ingest ──► index ──► lookup / search / ask
```

## Setup

```bash
Rscript -e 'install.packages(c("tidyverse","httr","stringr","readr","knitr","gridExtra","RColorBrewer"))'
pip install -r requirements.txt   # Python 3.10+
cp .env.example .env              # Windows: copy .env.example .env
```

Then open `.env` (`open -e .env` on Mac, `notepad .env` on Windows) and set
`GEMINI_API_KEY=your-key` and/or `OPENAI_API_KEY=…`. OpenAI is tried first.

`.env` is git-ignored and loaded automatically. If a setting appears twice, the **first** line
wins, so edit lines rather than appending new ones.

## Stage 1: R pipeline

**Input:** `MOLT4 mutations.csv`, the MOLT-4 mutation table from the
[DepMap portal](https://depmap.org/portal/) (GRCh38). It uses the columns `Gene`,
`Uniprot ID`, `Variant Info` and `Protein Change`.

```bash
Rscript MOLT4_protein_peptide_databases.R              # or add a path to the CSV
```

**Steps:**
1. Keep missense variants that have a UniProt ID (including `missense_variant&splice_region_variant`).
2. Fetch canonical sequences from UniProt. They're cached in `uniprot_sequence_cache.csv`.
3. Check the reference amino acid and apply the mutation. Variants that fail are logged with a reason.
4. Reverse each mutated protein to make a decoy.
5. Cut the peptide from the 2nd K/R upstream to the 2nd K/R downstream (or to the protein end).
6. Spot-check 5 random variants, then plot peptide lengths.

**Counts:**

| Input rows | With UniProt ID | Missense | In databases | Dropped | FASTA entries per database |
|---:|---:|---:|---:|---:|---:|
| 3826 | 2271 | 1871 | 1869 | 2 | 3738 (forward + decoy) |

The two dropped variants are HELZ2 R563L (UniProt has A, not R, at that position) and
TSC1 P1142_P1143delinsQT (not a single substitution).

**Outputs:**

| File | Contents |
|------|----------|
| `MOLT4_mutated_protein_database.fasta` | Variant proteins and their decoys |
| `MOLT4_mutated_peptide_database.fasta` | Variant peptides and their decoys |
| `MOLT4_peptide_data.csv` | Each peptide's sequence, start/end and mutation position |
| `MOLT4_mutated_protein_output.csv` | Canonical and mutated sequences plus DepMap columns |
| `MOLT4_mutations_with_sequences.csv` | Filtered variants plus canonical sequences |
| `MOLT4_dropped_variants.csv` | Dropped variants and the reason |
| `MOLT4_analysis_summary.csv` | The counts above |
| `MOLT4_variant_peptide_length_histogram.png` | Peptide length distribution |

FASTA header (decoys start with `Rev_sp`):
`>Fwd_spP750Q|Q5SV97-1|PERM1_P750Q OS=Homo sapiens GN=PERM1`

## Stage 2: Variant assistant

| Command | What it does | Needs |
|---------|--------------|-------|
| `python -m rag_assistant ingest` | Builds the index; downloads ClinVar (about 400 MB) once | Network |
| `python -m rag_assistant embed` | Adds embeddings for meaning-based search; resumes if stopped | API key |
| `python -m rag_assistant lookup NRAS G12C` | Exact lookup; also `--rsid rs…` or `--uniprot P…` | — |
| `python -m rag_assistant search "text"` | Keyword search | — |
| `python -m rag_assistant ask "question"` | LLM answer with sources | API key |
| `python -m rag_assistant ask "question" --no-llm` | The records `ask` would use, without the LLM | — |

Other options: `ingest --clinvar FILE` uses a local ClinVar file, `ask --provider gemini`
uses one provider only, and `ask --json` prints the full result.

**How `ask` works:**
1. An identifier in the question (gene + change such as `G12C` or `p.Gly12Cys`, an rsID, or a
   UniProt ID) is matched exactly. A gene on its own lists its variants. Anything else is
   searched.
2. Each variant comes with its linked records:
   - its peptide
   - its drop reason
   - its ClinVar match (by allele, then rsID, then protein change)
   - its mass spec hits, with `covers_variant` true only if the peptide spans the mutation.
3. The LLM answers from those records only, and can look up more. Every fact cites a source
   ID: `depmap:rowN`, `peptide:rowN`, `dropped:rowN`, `clinvar:VCV…` or `ms:rowN`.

**Example output:**

```
NRAS G12C is pathogenic/likely pathogenic in ClinVar [clinvar:VCV…] and a hotspot [depmap:row183].

[gemini gemini-3.6-flash, retrieval: exact]
```

- `retrieval: exact` means an identifier matched. `retrieval: search` means the answer is
  based on search results, so check it more carefully.
- `Warning: cited sources not in the retrieved records` means part of the answer is
  unsupported.
- For data that isn't there (for example BRAF V600E), the correct answer is "not found".

## Troubleshooting

| Message | Fix |
|---------|-----|
| `No API key for openai or gemini` | Create `.env` with your key (see Setup) |
| `404 … no longer available` | The model was retired. Set `GEMINI_CHAT_MODEL=` in `.env` to the model the error names |
| `503 … high demand` | The model is overloaded. Wait, or switch to another model |
| `APITimeoutError` | Retry, raise `LLM_TIMEOUT` (default 180 s), or switch model |
| `Function calling is not enabled` | Use a regular flash model, not `lite` or `gemma` |
| `Could not read ClinVar file … delete it` | Delete the file in `data/` and run `ingest` again |
| `Stopped early: … RateLimitError` | Run `embed` again; it resumes |
| A `.env` setting seems ignored | The key appears twice. Edit the first line |

Default models are `gemini-3.6-flash` / `gemini-embedding-001` (Gemini) and `gpt-4o-mini` /
`text-embedding-3-small` (OpenAI). To list the Gemini models your key can use:

```bash
KEY=$(grep '^GEMINI_API_KEY=' .env | cut -d= -f2-)
curl -sS -H "x-goog-api-key: $KEY" https://generativelanguage.googleapis.com/v1beta/models | grep '"name"'
```

## Data notes

- **Mass spec:** the hits come from a FragPipe search in
  [ndaw365/MOLT4-Checking](https://github.com/ndaw365/MOLT4-Checking). Only 7 of 41 hits
  span the mutated residue: EVL D20N, COPS7B A224T, TUBA4A E77D, CDC45 E259K, DUSP7 R102C,
  BRD4 P1131A and P2RX4 Y378F. The rest also match the normal protein. The search ran before
  59 `missense_variant&splice_region_variant` variants were added, so those 59 haven't been
  searched.
- **Histogram:** counts the 1765 forward peptides with at least 2 K/R on each side. 0.6% are
  under 7 aa, 76.4% are 7–50 aa and 22.9% are over 50 aa. The plot is cut off at 100 aa.
- **ClinVar:** "Vep Clin Sig" is DepMap's copy of ClinVar, while `clinvar:` entries come from
  the release you ingested. They can differ.

## Tests

```bash
python -m pytest          # 65 tests; fake ClinVar and LLM, no keys needed (CI runs them on every push)
python eval/run.py        # retrieval on 40 questions: 40/40 correct, ~0.3 ms each (CSV scan: 8/40)
python eval/run.py --llm  # also checks that LLM answers cite the right record (needs a key)
```

The evaluation questions come from the same data, so they test the lookup logic, not
open-ended questions.
