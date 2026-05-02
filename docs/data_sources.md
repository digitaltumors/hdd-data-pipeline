# Data Sources (HDD_v2)

This document lists the external inputs and generated datasets used to build the Harmonized Drug Dataset Version 2 (HDD_v2). Versions and URLs are defined in `config/pipeline.yaml`.

## External data sources

| Source | Version | URL / Endpoint | Access Method | Format | Notes |
| --- | --- | --- | --- | --- | --- |
| JUMP-CP cpg0016 MAE | v2.0.0 | Configured under `sub_dataset.jump.url` | Direct download (Snakemake `download_sub_dataset_mae`) | RDS | `metadata(mae)$Drug.Metadata` provides JUMP membership and source IDs. |
| OASIS cpg0037 MAE | v2.0.0 | Configured under `sub_dataset.oasis.url` | Direct download (Snakemake `download_sub_dataset_mae`) | RDS | `metadata(mae)$Drug.Metadata` provides OASIS membership and source IDs. |
| GEOM MAE | v2.0.0 | Configured under `sub_dataset.geom.url` | Direct download (Snakemake `download_sub_dataset_mae`) | RDS | `metadata(mae)$Drug.Metadata` provides GEOM membership and source SMILES. |
| LINCS MAE | v2.0.0 | Configured under `sub_dataset.lincs.url` | Direct download (Snakemake `download_sub_dataset_mae`) | RDS | `metadata(mae)$Drug.Metadata` provides LINCS membership and CMap names. |
| DeepChem BBBP | N/A | https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv | Direct download (Snakemake `download_DeepChem`) | CSV | Used to align compounds for annotation processing. |
| DeepChem ToxCast | N/A | https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/toxcast_data.csv.gz | Direct download (Snakemake `download_DeepChem`) | CSV (gzip) | Converted into experiment matrices. |
| DeepChem Tox21 | N/A | https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/tox21.csv.gz | Direct download (Snakemake `download_DeepChem`) | CSV (gzip) | Converted into experiment matrices. |
| DeepChem SIDER | N/A | https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/sider.csv.gz | Direct download (Snakemake `download_DeepChem`) | CSV (gzip) | Converted into experiment matrices. |
| DeepChem ClinTox | N/A | https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/clintox.csv.gz | Direct download (Snakemake `download_DeepChem`) | CSV (gzip) | Converted into experiment matrices. |
| BindingDB | 202512 | Configured from `binding_db.base_url`, `binding_db.subset`, and `binding_db.version` | Direct download (Snakemake `download_BindingDB`) | TSV (zip) | Filtered to human target records and converted into a target-by-compound affinity matrix. |
| AnnotationDB compound list | live | https://v2annotationdb.bhklab.ca/compound/all | API request (Snakemake `fetch_AnnotationDB_raw`) | JSON | Used to enumerate compounds and request detailed records. |
| AnnotationDB compound details | live | https://v2annotationdb.bhklab.ca/compound/many | Batched API requests (Snakemake `fetch_AnnotationDB_raw`) | JSON | Provides bioassays, mechanisms, and toxicity annotations. |

For license and citation requirements, consult each source website or associated publication.

## Generated datasets

| Dataset | Location | Created By | Inputs |
| --- | --- | --- | --- |
| AnnotationDB compact JSONL | `data/procdata/ANNOTATION_DB/compound_details.jsonl` | `workflow/scripts/fetch_annotationdb.py` | AnnotationDB API |
| Sub-dataset drug metadata | `data/procdata/sub_dataset/*_drug_metadata.tsv` | `workflow/scripts/extract_sub_dataset_drug_metadata.R` | Curated sub-dataset MAE RDS files |
| Sub-dataset parity report | `data/procdata/sub_dataset/parity/` | `workflow/scripts/process_annotationdb.py` | colData, extracted sub-dataset drug metadata |
| colData metadata | `data/procdata/colData.csv` | `workflow/scripts/process_annotationdb.py` | AnnotationDB JSONL, sub-dataset drug metadata, BBBP |
| Bioassay matrix | `data/procdata/experiments/bioassays.csv` | `workflow/scripts/process_annotationdb.py` | AnnotationDB JSONL |
| DeepChem experiment matrices | `data/procdata/experiments/{toxcast,tox21,sider,clintox}.csv` | `workflow/scripts/make_deepchem_experiments.py` | colData, DeepChem CSVs |
| BindingDB cleaned table | `data/procdata/BINDING_DB/BindingDB_*_cleaned.csv` | `workflow/scripts/process_bindingdb.py` | BindingDB bulk TSV zip |
| BindingDB affinity matrix | `data/procdata/experiments/binding_db.csv` | `workflow/scripts/make_bindingdb_experiments.py` | colData, cleaned BindingDB table |
| Morgan fingerprints | `data/procdata/experiments/fingerprints/Morgan.*.mtx` and `fingerprint_columns.tsv` | `workflow/scripts/make_fingerprints.py` | parseable colData SMILES |
| HDD_v2 MAE | `data/results/HDD_v2.RDS` | `workflow/scripts/construct_MAE.R` | colData + experiment matrices |
| HDD_v2 CSV exports | `data/results/HDD_v2_csv/` | `workflow/scripts/export_mae_csvs.R` | HDD_v2 MAE |
