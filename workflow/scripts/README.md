# Workflow Scripts

This directory contains the executable scripts used by the Snakemake workflow to build HDD_v2.3.

## Script catalog

- `fetch_annotationdb.py`
  - Fetches compound metadata from AnnotationDB (`/compound/all` and `/compound/many`) and writes a compact JSONL intermediate.
  - Output: `data/procdata/ANNOTATION_DB/compound_details.jsonl`.

- `process_annotationdb.py`
  - Parses the AnnotationDB JSONL, flattens source-specific toxicity metadata, and joins extracted sub-dataset drug metadata plus BBBP metadata.
  - Outputs: `data/procdata/colData.csv` and `data/procdata/experiments/bioassays.csv`.

- `extract_sub_dataset_drug_metadata.R`
  - Reads drug metadata from downloaded curated MAE and PharmacoSet RDS objects and normalizes their compound identity fields.
  - Output: `data/procdata/sub_dataset/*_drug_metadata.tsv`.

- `make_deepchem_experiments.py`
  - Reshapes DeepChem task datasets into `HDD.Compound.ID`-by-assay matrices.
  - Outputs: `data/procdata/experiments/{toxcast,tox21,sider,clintox}.csv`.

- `make_fingerprints.py`
  - Generates Morgan count fingerprints from parseable SMILES for configured radii and dimensions.
  - Output: `data/procdata/experiments/fingerprints/Morgan.*.mtx` and `fingerprint_columns.tsv`.

- `construct_MAE.R`
  - Assembles all experiment matrices and colData into a `MultiAssayExperiment`.
  - Output: `data/results/HDD_v2.3.RDS`.

- `export_mae_csvs.R`
  - Exports MAE-backed CSVs for parity with the RDS output, including sparse fingerprint assays.
  - Output: `data/results/HDD_v2.3_csv/`.

- `archive_mae_csvs.py`
  - Archives `data/results/HDD_v2.3_csv/` into the `tar.gz`.
  - Output: `data/results/HDD_v2.3_csv.tar.gz`.

## Notes

- Scripts are invoked by rules in `workflow/rules/` and use paths from `damply.dirs`.
- If you change inputs or URLs in `config/pipeline.yaml`, re-run the pipeline to regenerate outputs.
