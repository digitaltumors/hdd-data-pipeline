# Usage Guide

## Project configuration

All pipeline settings live in `config/pipeline.yaml`. The most common edits are:

- **Data source versions and URLs**
  - BindingDB version and subset
  - LINCS compound info release
  - JUMP-CP metadata release
  - DeepChem dataset URLs
  - AnnotationDB endpoint
- **BindingDB filtering**
  - Organism allowlist
  - Columns retained in the cleaned export
- **Fingerprint parameters**
  - Morgan radii and vector dimensions

If you change versions or URLs, update `docs/data_sources.md` so the provenance stays current.

## Data locations

The pipeline writes data into three main locations:

- `data/rawdata/`: raw downloads (BindingDB archive/TSV, DeepChem tables, LINCS, JUMP-CP).
- `data/procdata/`: processed datasets (AnnotationDB compact JSONL, BindingDB cleaned table, colData, experiments, sparse fingerprints).
- `data/results/`: final HDD_v1 output (`HDD_v1.RDS`).
- `data/results/HDD_v1_csv/`: MAE-derived CSVs for `colData` and dense assays. Sparse fingerprint assays remain in `data/procdata/experiments/fingerprints/` as `.mtx`.

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

The pipeline also writes MAE-derived CSV exports to `data/results/HDD_v1_csv/`.

To bundle those exports with the sparse fingerprint `.mtx` assays, run:

```bash
pixi run zip-output
```

This writes `data/results/HDD_v1_csv.tar.gz`.

## Quality control

Render the QC report after the pipeline has produced `HDD_v1.RDS`:

```bash
pixi run knit-qc
```

The report is saved to `qc/hdd_quality_control.html`.
