# Developer Notes

## Scope Decisions

- HDD v2 is a compound-level `MultiAssayExperiment` assembled from AnnotationDB plus curated sub-dataset MAE inputs.
- JUMP-CP, OASIS, GEOM, and LINCS membership is sourced from locked curated MAE RDS files configured under `sub_dataset`.
- The pipeline reads `metadata(mae)$Drug.Metadata` from each sub-dataset MAE and uses those tables to populate `In.JUMP.CP`, `In.OASIS`, `In.GEOM`, `In.LINCS`, and their source key columns.
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
- `SMILES` is the HDD working structure column. It stores the best available structure string for each HDD row.
- `AnnotationDB.SMILES` is the structure string specifically sourced from AnnotationDB. It is populated only for compounds that map to AnnotationDB.
- `GEOM.Source.SMILES` is the GEOM source compound key retained for GEOM membership, traceability, and parity checks. It is not the general HDD structure column.
- For GEOM-derived rows, `SMILES` is filled from the first available value among `GEOM.Canonical.SMILES`, `GEOM.RDKit.SMILES`, `GEOM.Source.SMILES`, and `AnnotationDB.SMILES` when the row does not already have an AnnotationDB structure.
- AnnotationDB-only fields remain missing for source-only compounds that do not map to AnnotationDB.

## Assay Decisions

- DeepChem tasks (ToxCast, Tox21, SIDER, ClinTox) are converted into `HDD.Compound.ID`-by-assay matrices for consistent MAE ingestion.
- BindingDB is included as a target-by-compound assay keyed by `HDD.Compound.ID`.
- BindingDB records are matched to HDD compounds by `Pubchem.CID`, so HDD compounds without PubChem CIDs are retained in `colData` but absent from the BindingDB assay.
- The BindingDB value is the minimum exact numeric `Ki (nM)` or `Kd (nM)` reported for each human target-compound pair after filtering out PubChem BioAssay-linked records.
- Morgan count fingerprints use the HDD `SMILES` column as input to RDKit.
- Morgan count fingerprints are generated only for compounds with parseable `SMILES` and stored as sparse Matrix Market experiments with `fingerprint_columns.tsv` mapping matrix columns back to `HDD.Compound.ID`.
- Compounds with no assay columns and no parseable `SMILES` remain in MAE `colData` for metadata completeness, even though they are absent from `sampleMap` until represented in an assay.

## Output Decisions

- The primary pipeline output is `data/results/HDD_v2.RDS`.
- MAE-derived CSV exports are written under `data/results/HDD_v2_csv/`.
- Sparse fingerprint assays stay in Matrix Market format in the table archive, with `fingerprint_columns.tsv` included beside the `.mtx` files.
