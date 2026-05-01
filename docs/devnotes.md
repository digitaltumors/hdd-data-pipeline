# Developer Notes

## 2026-04-30 - Dataset naming for HDD_v2

- The pipeline output is now `HDD_v2.RDS`, with CSV exports under `data/results/HDD_v2_csv/`.
- LINCS membership is represented as `In.LINCS` and `LINCS.ID`; previous release-specific names are not emitted.
- OASIS and GEOM membership are joined from published minimal membership artifacts by exact InChIKey.
- Legacy affinity-source rules and scripts were removed from the public HDD v2 build.

## 2025-12-24 - Previous release naming

- The output file remains a `MultiAssayExperiment` assembled by `workflow/scripts/construct_MAE.R`.

## 2025-12-24 - Compound universe and metadata strategy

- The compound universe is sourced from AnnotationDB (`/compound/all`), then enriched with detailed records from `/compound/many`.
- LINCS and JUMP-CP compound metadata are merged to harmonize identifiers used across public resources.
- DeepChem BBBP data is used during metadata processing to align SMILES and CID mappings.

## 2025-12-24 - Assay integration decisions

- DeepChem tasks (ToxCast, Tox21, SIDER, ClinTox) are converted into CID-by-assay matrices for consistent MAE ingestion.
- Morgan count fingerprints (configurable radii and dimensions) are generated from colData SMILES and stored as sparse Matrix Market experiments.
