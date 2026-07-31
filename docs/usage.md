# Usage Guide

## Project configuration

All pipeline settings live in `config/pipeline.yaml`. The most common edits are:

- **Data source versions and URLs**
  - Curated JUMP-CP, OASIS, GEOM, and LINCS MAE release URLs
  - CTRPv2 and NCI60-2026 PharmacoSet release URLs
  - DeepChem dataset URLs
  - AnnotationDB endpoint
- **Fingerprint parameters**
  - Morgan radii and vector dimensions

If you change versions or URLs, update `docs/data_sources.md` so the provenance stays current.

ATC data from AnnotationDB requires an authorization key. Copy `.env.example`
to `.env` and set `ANNOTATIONDB_API_KEY`, or set the same variable in the
process environment. The `.env` file is ignored by Git. If the key is absent,
the pipeline still runs and fetches the other AnnotationDB fields, but excludes
ATC data.

## Data locations

The pipeline writes data into three main locations:

- `data/rawdata/`: raw downloads (DeepChem tables and curated MAE and PharmacoSet objects).
- `data/procdata/`: processed datasets (AnnotationDB compact JSONL, extracted sub-dataset metadata, colData, experiments, sparse fingerprints, and the fingerprint column map).
- `data/results/`: final HDD_v2.3 output (`HDD_v2.3.RDS`).
- `data/results/HDD_v2.3_csv/`: MAE-derived CSVs for `colData` and all assay matrices, including Morgan fingerprints.

Raw and processed files are not tracked in Git, so make sure you archive them externally if you need to preserve a run.

## Running the pipeline

Install dependencies with Pixi:

```bash
pixi install
```

Run Snakemake from the repository root:

```bash
pixi run snakemake -c 1
```

The pipeline also writes MAE-derived CSV exports to `data/results/HDD_v2.3_csv/`.

To bundle those exports, run:

```bash
pixi run zip-output
```

This writes `data/results/HDD_v2.3_csv.tar.gz`.

## Inspecting the RDS

Load `MultiAssayExperiment` before inspecting the RDS so S4 generics are resolved
through the attached package rather than namespace-qualified one-liners:

```r
suppressPackageStartupMessages(library(MultiAssayExperiment))

hdd_mae <- readRDS("data/results/HDD_v2.3.RDS")
validObject(hdd_mae)
dim(colData(hdd_mae))
names(experiments(hdd_mae))
```

## Quality control

Render the QC report after the pipeline has produced `HDD_v2.3.RDS`:

```bash
pixi run knit-qc
```

The report is saved to `qc/hdd_quality_control.html`.
