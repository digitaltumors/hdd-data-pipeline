# Developer Notes

## Scope Decisions

- HDD is a compound-level `MultiAssayExperiment` assembled from AnnotationDB plus curated sub-dataset inputs.
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
- DeepChem BBBP is joined during metadata processing by exact SMILES to fill `BBB.Permeable`. NCATS Tox21 uses a unique full InChIKey match against the assembled HDD compound universe. Neither needs an external SMILES-to-CID mapping.
- ChEMBL drug indications from AnnotationDB are retained as a one-to-many table in `metadata(mae)$Drug.Indications` rather than flattened into `colData`.
- AnnotationDB toxicity metadata is flattened into `colData` as LTKB label-set fields, LiverTox fields, DIRIL fields, and DICT fields when returned by the API. The older collapsed DILI/hepatotoxicity columns are not emitted.
- AnnotationDB `atc_code` is flattened into `colData` as `ATC.Code`. The field is included when `ANNOTATIONDB_API_KEY` is set in `.env`; otherwise the pipeline continues and `ATC.Code` remains missing.
- `SMILES` is the HDD working structure column. It stores the best available structure string for each HDD row.
- `AnnotationDB.SMILES` is the structure string specifically sourced from AnnotationDB. It is populated only for compounds that map to AnnotationDB.
- `GEOM.Source.SMILES` is the GEOM source compound key retained for GEOM membership, traceability, and parity checks. It is not the general HDD structure column.
- For GEOM-derived rows, `SMILES` is filled from the first available value among `GEOM.Canonical.SMILES`, `GEOM.RDKit.SMILES`, `GEOM.Source.SMILES`, and `AnnotationDB.SMILES` when the row does not already have an AnnotationDB structure.
- AnnotationDB-only fields remain missing for source-only compounds that do not map to AnnotationDB.

## Assay Decisions

- The workflow no longer uses the DeepChem SIDER and Tox21 inputs. It omits the SIDER experiment and re-curates Tox21 from the original NCATS 2014 Challenge training archive. DeepChem BBBP remains a compound metadata input. ClinTox and ToxCast remain deprecated.
- Tox21 decisions are per HDD compound and assay. Only explicit `0` (inactive) and `1` (active) are votes. If the non-missing calls agree, keep that label; if both occur, set that cell to `NA`. No call also yields `NA`. A conflict does not remove the compound or its other assay labels.
- Do not pick a first row, let any positive win, or majority-vote: sample order is arbitrary and the binary calls do not say which sample is more reliable.
- Require a unique full InChIKey match to HDD. Leave unmatched and ambiguous samples unassigned and report them. Historical `DSSTox_CID` is provenance, not an HDD join key.
- Preserve the `Tox21` experiment and `Tox21.tsv` assay-by-compound interface.
- Morgan count fingerprints use the HDD `SMILES` column as input to RDKit.
- Morgan count fingerprints are generated only for compounds with parseable `SMILES` and stored in the MAE as sparse assays.
- Compounds with no assay columns and no parseable `SMILES` remain in MAE `colData` for metadata completeness, even though they are absent from `sampleMap` until represented in an assay.

## Output Decisions

- The primary pipeline output is `data/results/<dataset_id>.RDS`, where `<dataset_id>` is `dataset.id` in `config/pipeline.yaml`.
- MAE-derived flat exports are written under `data/results/<dataset_id>_tables/`.
- Dense tables are TSV; sparse Morgan fingerprint assays use Matrix Market files with explicit feature and compound mappings.
- `data/results/<dataset_id>_tables.tar.gz` is the shareable gzip-compressed archive.

## DailyMed mapping

DailyMed travels through the existing AnnotationDB JSONL, compound processor
and MAE constructor. No separate download or prior-release enrichment branch
is needed. Both the response CID and label PubChem CID must match exactly;
there is no salt/ingredient propagation or second label selection. Conflicting
duplicate responses and unexpected label cardinality fail rather than silently
selecting a record.

`DailyMed.Query.Status` distinguishes a returned label, a returned compound
without a label, a CID absent from AnnotationDB responses, and no PubChem CID.
The compact cache must contain the label key even when its value is null;
older caches without it require regeneration, not reinterpretation as missing
labels. The intermediate TSV uses `\N` only for DailyMed nulls, preserving
empty strings and multiline text. The R constructor reads those columns with
their own null convention; existing compound-column parsing stays unchanged.

Keep one QC source and HTML, update documentation in place, and remove retired
workflow branches. Developer checks belong outside the published report;
historical source snapshots retain their retrieval metadata without dates in
filenames. Existing curated releases are preserved when changing the workflow.
