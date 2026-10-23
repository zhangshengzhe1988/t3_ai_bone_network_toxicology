# AI-associated bone loss: network toxicology of aromatase inhibitors

Code and derived data for the manuscript:

> **Aromatase inhibitors as endocrine-disrupting pharmaceuticals: network toxicology links CYP19A1 inhibition to human bone risk and cross-species effects**
> Shengzhe Zhang, Han Xu, Ting Gu, Hua Dai (Yangzhou University)

## What this repository contains

| Folder | Content |
|--------|---------|
| `code/` | Analysis and figure scripts (Python 3.10 + R 4.6 for WGCNA/limma) |
| `results/tables/` | Derived summary tables reported in the manuscript and Supplementary Tables |
| `results/figures/` | Main figures (Fig. 1–7) and supplementary figures (S1–S11), vector PDF |
| `data/` | Data-source notes (raw data are public and not redistributed here) |
| `manuscript/` | Manuscript status; text will be deposited after acceptance |

## Reproducing the analysis

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# 1. FAERS disproportionality (openFDA API)
python code/t3_manuscript_v2.py        # builds manuscript and figure legends

# 2. Figures
python code/assemble_figures.py        # assembles Fig. 1–7 (A4, 210 mm)

# 3. Exposure–activity ratio
python code/t3_ear_exposure.py
```

WGCNA and limma are run in R; see `code/S7_wgcna.R` and the per-cohort limma scripts.
Seeds, thresholds and software versions are listed in the manuscript Methods (section 2).

## Data sources

| Source | Access |
|--------|--------|
| FAERS / openFDA | https://open.fda.gov/apis/drug/event/ (queried 2026-10-13, window 2013-01-01 to 2026-06-30) |
| GEO GSE35958, GSE13850, GSE56814, GSE56815, GSE201774 | https://www.ncbi.nlm.nih.gov/geo/ |
| ChEMBL | https://www.ebi.ac.uk/chembl/ |
| SwissTargetPrediction | http://www.swisstargetprediction.ch/ |
| GeneCards | https://www.genecards.org/ |
| STRING v12 | https://string-db.org/ (score ≥ 0.4) |
| PDB 3S7S, AlphaFold2 models | https://www.rcsb.org/ , https://alphafold.ebi.ac.uk/ |
| PubChem | https://pubchem.ncbi.nlm.nih.gov/ |

No raw data are redistributed; all sources are public and versioned by their accessions.

## Licence

Code: MIT (`LICENSE`). Derived tables and figures: CC BY 4.0.

## Citation

See `CITATION.cff`. A Zenodo DOI will be added on first release.
