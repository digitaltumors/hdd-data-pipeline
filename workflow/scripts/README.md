# Workflow Scripts

This directory contains executable scripts used by Snakemake to build HDD v3.

## Script catalog

- `fetch_annotationdb.py`
  - Fetches compound metadata from AnnotationDB (`/compound/all` and `/compound/many`) and writes a compact JSONL intermediate.
  - Output: `data/procdata/ANNOTATION_DB/compound_details.jsonl`.

- `process_annotationdb.py`
  - Parses the AnnotationDB JSONL, flattens source-specific toxicity metadata, and joins extracted sub-dataset drug metadata plus BBBP metadata.
  - Outputs: `colData.tsv`, `bioassays.tsv`, `drug_indications.tsv`, and source-key parity reports.

- `extract_sub_dataset_drug_metadata.R`
  - Reads drug metadata from downloaded curated MAE and PharmacoSet RDS objects and normalizes their compound identity fields.
  - Output: `data/procdata/sub_dataset/*_drug_metadata.tsv`.

- `make_deepchem_experiments.py`
  - Matches DeepChem structures by exact SMILES and local full InChIKey,
    applies the AnnotationDB preference and first-source-row collision policy,
    and reshapes accepted records into `HDD.Compound.ID`-by-assay matrices.
  - Outputs: `data/procdata/experiments/{tox21,sider}.tsv` and DeepChem match
    audits under `data/procdata/metadata/`.

- `make_fingerprints.py`
  - Generates Morgan count fingerprints from parseable SMILES for configured radii and dimensions.
  - Output: `data/procdata/experiments/fingerprints/Morgan.*.mtx` and `fingerprint_columns.tsv`.

- `construct_MAE.R`
  - Assembles all experiment matrices and colData into a `MultiAssayExperiment`.
  - Output: `data/results/HDD_v3.RDS`.

- `export_mae_tables.R`
  - Exports dense TSVs and sparse Matrix Market assays without densification.
  - Output: `data/results/HDD_v3_tables/`.

- `validate_release.R`
  - Enforces object identity, exact experiment allowlist, DeepChem dimensions,
    deprecated-assay absence, indication count and foreign keys, and
    `sampleMap` integrity.

## Notes

- Scripts are invoked by rules in `workflow/rules/`; all paths and controls come from `config/pipeline.yaml`.
- If you change inputs or URLs in `config/pipeline.yaml`, re-run the pipeline to regenerate outputs.
