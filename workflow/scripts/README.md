# Workflow Scripts

This directory contains the executable scripts used by the Snakemake workflow to build HDD_v1.

## Script catalog

- `fetch_annotationdb.py`
  - Downloads compound metadata from AnnotationDB (`/compound/all` and `/compound/many`).
  - Output: `data/rawdata/ANNOTATION_DB/compound_details.jsonl`.

- `process_annotationdb.py`
  - Parses the AnnotationDB JSONL and joins LINCS, JUMP-CP, and BBBP metadata.
  - Outputs: `data/procdata/colData.csv` and `data/procdata/experiments/bioassays.csv`.

- `make_bindingdbd_experiments.py`
  - Converts the cleaned BindingDB table into a CID-by-target affinity matrix.
  - Output: `data/procdata/experiments/binding_db.csv`.

- `make_deepchem_experiments.py`
  - Reshapes DeepChem task datasets into CID-by-assay matrices.
  - Outputs: `data/procdata/experiments/{toxcast,tox21,sider,clintox}.csv`.

- `make_fingerprints.py`
  - Generates Morgan count fingerprints from SMILES for configured radii and dimensions.
  - Output: `data/procdata/experiments/fingerprints/Morgan.*.csv`.

- `construct_MAE.R`
  - Assembles all experiment matrices and colData into a `MultiAssayExperiment`.
  - Output: `data/results/HDD_v1.RDS`.

- `export_mae_csvs.R`
  - Exports MAE-backed CSVs for parity with the RDS output.
  - Output: `data/results/HDD_v1_csv/`.

- `utils.py`
  - Helper functions shared by AnnotationDB processing scripts (assay filters, field parsing, etc.).

- `make_colData.py`
  - Legacy helper for colData generation. Not used by the current Snakemake workflow.

## Notes

- Scripts are invoked by rules in `workflow/rules/` and use paths from `damply.dirs`.
- If you change inputs or URLs in `config/pipeline.yaml`, re-run the pipeline to regenerate outputs.
