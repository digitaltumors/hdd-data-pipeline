# Developer Notes

## 2026-05-01 - Sub-dataset MAE membership inputs

- JUMP-CP, OASIS, GEOM, and LINCS membership are sourced from locked curated MAE RDS files configured under `sub_dataset`.
- The pipeline reads `metadata(mae)$Drug.Metadata` from each sub-dataset MAE and uses those tables to populate `In.JUMP.CP`, `In.OASIS`, `In.GEOM`, `In.LINCS`, and their source key columns.
- `Pubchem.CID` is nullable and is no longer the primary HDD key. `HDD.Compound.ID` is generated from an internal SHA-256 identity hash and is used as the MAE colData row name and assay column key.
- Non-AnnotationDB compounds from sub-dataset metadata are retained as normal HDD rows with missing AnnotationDB-only fields.
- `OASIS.ID` is kept as reported by the OASIS MAE. OASIS rows with missing `OASIS.ID` are retained and reported separately in the parity report.
- `LINCS.CMap.Name` is the LINCS source key. `LINCS.ID` is not emitted.

## 2026-04-30 - Dataset naming for HDD_v2

- The pipeline output is now `HDD_v2.RDS`, with CSV exports under `data/results/HDD_v2_csv/`.
- Legacy affinity-source rules and scripts were removed from the public HDD v2 build.

## 2025-12-24 - Previous release naming

- The output file remains a `MultiAssayExperiment` assembled by `workflow/scripts/construct_MAE.R`.

## 2025-12-24 - Compound universe and metadata strategy

- The compound universe is sourced from AnnotationDB (`/compound/all`) plus curated sub-dataset MAE drug metadata, then enriched with detailed AnnotationDB records from `/compound/many` where available.
- DeepChem BBBP data is used during metadata processing to align SMILES and CID mappings.

## 2025-12-24 - Assay integration decisions

- DeepChem tasks (ToxCast, Tox21, SIDER, ClinTox) are converted into `HDD.Compound.ID`-by-assay matrices for consistent MAE ingestion.
- Morgan count fingerprints (configurable radii and dimensions) are generated from colData SMILES and stored as sparse Matrix Market experiments.
