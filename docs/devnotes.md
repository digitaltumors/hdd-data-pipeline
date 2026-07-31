# Developer Notes

## Scope Decisions

- HDD v2.3 is a compound-level `MultiAssayExperiment` assembled from AnnotationDB plus curated sub-dataset inputs.
- JUMP-CP, OASIS, GEOM, and LINCS membership is sourced from locked curated MAE RDS files configured under `sub_dataset`.
- CTRPv2 and NCI60 membership is sourced from downloaded PharmacoSet RDS files configured under `sub_dataset`.
- The pipeline reads `metadata(mae)$Drug.Metadata` from each sub-dataset MAE and the treatment metadata from each PharmacoSet to populate `In.JUMP.CP`, `In.OASIS`, `In.GEOM`, `In.LINCS`, `In.CTRP`, `In.NCI60`, and their source key columns.
- Every sub-dataset input is downloaded from a configured URL; the release pipeline does not depend on sibling-repository or machine-local paths.
- Non-AnnotationDB compounds from sub-dataset metadata are retained as normal HDD rows with missing AnnotationDB-only fields.
- Source-specific assay rules are kept only when they produce MAE assays keyed by `HDD.Compound.ID`.

## Identifier Decisions

- `Pubchem.CID` is nullable and is not the primary HDD key.
- `HDD.Compound.ID` is generated from an internal SHA-256 identity hash and is used as the MAE `colData` row name and assay column key.
- The hash identity priority is PubChem CID, then InChIKey, then source-specific identity. The identity string is used only internally and is not emitted as a public column.
- `OASIS.ID` is kept as reported by the OASIS MAE. OASIS rows with missing `OASIS.ID` are retained and reported separately in the parity report.
- `LINCS.CMap.Name` is the LINCS source key. `LINCS.ID` is not emitted.
- Public columns with HDD-shared names are intended to be directly joinable across HDD-facing datasets; source-specific fields use source-specific prefixes.

## Metadata Decisions

- The compound universe is sourced from AnnotationDB (`/compound/all`) plus curated sub-dataset MAE drug metadata, then enriched with detailed AnnotationDB records from `/compound/many` where available.
- DeepChem BBBP data is used during metadata processing to align SMILES and CID mappings.
- AnnotationDB toxicity metadata is flattened into `colData` as LTKB label-set fields, LiverTox fields, DIRIL fields, and DICT fields when returned by the API. The older collapsed DILI/hepatotoxicity columns are not emitted.
- AnnotationDB `atc_code` is flattened into `colData` as `ATC.Code`. The field is included when `ANNOTATIONDB_API_KEY` is set in `.env` or the process environment; otherwise the pipeline continues and `ATC.Code` remains missing.
- `SMILES` is the HDD working structure column. It stores the best available structure string for each HDD row.
- `AnnotationDB.SMILES` is the structure string specifically sourced from AnnotationDB. It is populated only for compounds that map to AnnotationDB.
- `GEOM.Source.SMILES` is the GEOM source compound key retained for GEOM membership, traceability, and parity checks. It is not the general HDD structure column.
- For GEOM-derived rows, `SMILES` is filled from the first available value among `GEOM.Canonical.SMILES`, `GEOM.RDKit.SMILES`, `GEOM.Source.SMILES`, and `AnnotationDB.SMILES` when the row does not already have an AnnotationDB structure.
- AnnotationDB-only fields remain missing for source-only compounds that do not map to AnnotationDB.

## Assay Decisions

- DeepChem tasks (ToxCast, Tox21, SIDER, ClinTox) are converted into `HDD.Compound.ID`-by-assay matrices for consistent MAE ingestion.
- Morgan count fingerprints use the HDD `SMILES` column as input to RDKit.
- Morgan count fingerprints are generated only for compounds with parseable `SMILES` and stored in the MAE as sparse assays.
- Compounds with no assay columns and no parseable `SMILES` remain in MAE `colData` for metadata completeness, even though they are absent from `sampleMap` until represented in an assay.

## Output Decisions

- The primary pipeline output is `data/results/HDD_v2.3.RDS`.
- MAE-derived CSV exports are written under `data/results/HDD_v2.3_csv/`.
- Morgan fingerprint assays are exported as assay CSVs alongside the dense assay exports.
