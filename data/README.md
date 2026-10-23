# Data sources

Raw data are public and are **not** redistributed in this repository.

- **FAERS**: openFDA drug/event endpoint, queried 2026-10-13, window 2013-01-01 to 2026-06-30.
  The exact API calls are in `code/t3_manuscript_v2.py` (Methods 2.1).
- **GEO**: GSE35958, GSE13850, GSE56814, GSE56815, GSE201774.
- **Targets**: ChEMBL (merged), SwissTargetPrediction (p >= 0.1), GeneCards (osteoporosis query).
- **PPI**: STRING v12, combined score >= 0.4.
- **Structures**: PDB 3S7S (human CYP19A1); AlphaFold2 models AF-O42145-F1 (zebrafish Cyp19a1a),
  AF-Q96QK1-F1 (human VPS35).
- **Physicochemistry**: PubChem logKow; prescribing information (DailyMed) for therapeutic concentrations.

`results/tables/` holds derived summaries only; per-sample expression matrices are omitted
because they can be regenerated from the GEO accessions.
