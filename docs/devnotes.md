# Developer Notes

## 2025-12-24 - Dataset naming for HDD_v1

- The pipeline output was renamed to `HDD_v1.RDS` and all documentation now refers to the Harmonized Drug Dataset Version 1 (HDD_v1).
- Added an export step that regenerates MAE-derived CSVs under `data/results/HDD_v1_csvs/` for parity with the RDS output.
- The output file remains a `MultiAssayExperiment` assembled by `workflow/scripts/construct_MAE.R`.

## 2025-12-24 - Compound universe and metadata strategy

- The compound universe is sourced from AnnotationDB (`/compound/all`), then enriched with detailed records from `/compound/many`.
- LINCS and JUMP-CP compound metadata are merged to harmonize identifiers used across public resources.
- DeepChem BBBP data is used during metadata processing to align SMILES and CID mappings.

## 2025-12-24 - Assay integration decisions

- BindingDB is deprecated for now due to scope changes driven by data quality concerns.
- BindingDB is filtered to human targets and assays without PubChem AIDs to avoid mixing external bioassays with AnnotationDB AIDs.
- DeepChem tasks (ToxCast, Tox21, SIDER, ClinTox) are converted into CID-by-assay matrices for consistent MAE ingestion.
- Morgan count fingerprints (configurable radii and dimensions) are generated from colData SMILES and stored as separate experiments.
