# Data Sources (HDD_v1.1)

This document lists the external inputs and generated datasets used to build the Harmonized Drug Dataset Version 1.1 (HDD_v1.1). Versions and URLs are defined in `config/pipeline.yaml`.

## External data sources

| Source | Version | URL / Endpoint | Access Method | Format | Notes |
| --- | --- | --- | --- | --- | --- |
| BindingDB (All subset) | 202512 | https://www.bindingdb.org/rwd/bind/downloads/BindingDB_All_202512_tsv.zip | Direct download (Snakemake `download_BindingDB`) | TSV inside ZIP | Filtered to human targets and assays without PubChem AIDs. License and citation are per BindingDB. |
| LINCS compound info | 2020 | https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/compoundinfo_beta.txt | Direct download (Snakemake `download_LINCS`) | Tab-delimited text | Used for compound metadata cross-references. |
| JUMP-CP compound metadata | cpg0016 | https://github.com/jump-cellpainting/datasets/raw/refs/heads/main/metadata/compound.csv.gz | Direct download (Snakemake `download_JUMPCP`) | CSV (gzip) | Used for compound metadata alignment. |
| DeepChem BBBP | N/A | https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv | Direct download (Snakemake `download_DeepChem`) | CSV | Used to align compounds for annotation processing. |
| DeepChem ToxCast | N/A | https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/toxcast_data.csv.gz | Direct download (Snakemake `download_DeepChem`) | CSV (gzip) | Converted into experiment matrices. |
| DeepChem Tox21 | N/A | https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/tox21.csv.gz | Direct download (Snakemake `download_DeepChem`) | CSV (gzip) | Converted into experiment matrices. |
| DeepChem SIDER | N/A | https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/sider.csv.gz | Direct download (Snakemake `download_DeepChem`) | CSV (gzip) | Converted into experiment matrices. |
| DeepChem ClinTox | N/A | https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/clintox.csv.gz | Direct download (Snakemake `download_DeepChem`) | CSV (gzip) | Converted into experiment matrices. |
| AnnotationDB compound list | live | https://annotationdb.bhklab.ca/compound/all | API request (Snakemake `fetch_AnnotationDB_raw`) | JSON | Used to enumerate compounds and request detailed records. |
| AnnotationDB compound details | live | https://annotationdb.bhklab.ca/compound/many | Batched API requests (Snakemake `fetch_AnnotationDB_raw`) | JSON | Provides bioassays, mechanisms, and toxicity annotations. |

For license and citation requirements, consult each source website or associated publication.

## Generated datasets

| Dataset | Location | Created By | Inputs |
| --- | --- | --- | --- |
| BindingDB cleaned table | `data/procdata/BINDING_DB/BindingDB_All_202512_cleaned.csv` | `workflow/rules/processBindingDB.smk` | BindingDB ZIP archive |
| AnnotationDB compact JSONL | `data/procdata/ANNOTATION_DB/compound_details.jsonl` | `workflow/scripts/fetch_annotationdb.py` | AnnotationDB API |
| colData metadata | `data/procdata/colData.csv` | `workflow/scripts/process_annotationdb.py` | AnnotationDB JSONL, LINCS, JUMP-CP, BBBP |
| Bioassay matrix | `data/procdata/experiments/bioassays.csv` | `workflow/scripts/process_annotationdb.py` | AnnotationDB JSONL |
| BindingDB experiment matrix | `data/procdata/experiments/binding_db.csv` | `workflow/scripts/make_bindingdbd_experiments.py` | colData, BindingDB cleaned table |
| DeepChem experiment matrices | `data/procdata/experiments/{toxcast,tox21,sider,clintox}.csv` | `workflow/scripts/make_deepchem_experiments.py` | colData, DeepChem CSVs |
| Morgan fingerprints | `data/procdata/experiments/fingerprints/Morgan.*.mtx` | `workflow/scripts/make_fingerprints.py` | colData SMILES |
| HDD_v1.1 MAE | `data/results/HDD_v1.1.RDS` | `workflow/scripts/construct_MAE.R` | colData + experiment matrices |
| HDD_v1.1 CSV exports | `data/results/HDD_v1.1_csv/` | `workflow/scripts/export_mae_csvs.R` | HDD_v1.1 MAE |
