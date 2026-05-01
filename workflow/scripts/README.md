# Workflow Scripts

This directory contains the executable scripts used by the Snakemake workflow to build HDD_v2.

## Script catalog

- `fetch_annotationdb.py`
  - Fetches compound metadata from AnnotationDB (`/compound/all` and `/compound/many`) and writes a compact JSONL intermediate.
  - Output: `data/procdata/ANNOTATION_DB/compound_details.jsonl`.

- `process_annotationdb.py`
  - Parses the AnnotationDB JSONL and joins extracted sub-dataset drug metadata plus BBBP metadata.
  - Outputs: `data/procdata/colData.csv` and `data/procdata/experiments/bioassays.csv`.

- `extract_sub_dataset_drug_metadata.R`
  - Reads `metadata(mae)$Drug.Metadata` from each curated sub-dataset MAE RDS.
  - Output: `data/procdata/sub_dataset/*_drug_metadata.tsv`.

- `make_deepchem_experiments.py`
  - Reshapes DeepChem task datasets into `HDD.Compound.ID`-by-assay matrices.
  - Outputs: `data/procdata/experiments/{toxcast,tox21,sider,clintox}.csv`.

- `make_fingerprints.py`
  - Generates Morgan count fingerprints from SMILES for configured radii and dimensions.
  - Output: `data/procdata/experiments/fingerprints/Morgan.*.mtx`.

- `construct_MAE.R`
  - Assembles all experiment matrices and colData into a `MultiAssayExperiment`.
  - Output: `data/results/HDD_v2.RDS`.

- `export_mae_csvs.R`
  - Exports MAE-backed CSVs for parity with the RDS output, skipping sparse fingerprint assays.
  - Output: `data/results/HDD_v2_csv/`.

- `archive_mae_csvs.py`
  - Archives `data/results/HDD_v2_csv/` and injects sparse fingerprint `.mtx` assays into the `tar.gz` under `HDD_v2_csv/assays/`.
  - Output: `data/results/HDD_v2_csv.tar.gz`.

## Notes

- Scripts are invoked by rules in `workflow/rules/` and use paths from `damply.dirs`.
- If you change inputs or URLs in `config/pipeline.yaml`, re-run the pipeline to regenerate outputs.
