# HDD_v3 Data Pipeline

**Authors:** [James Bannon](https://github.com/jbannon), Michael Tran, Matthew Boccalon, Sisira Kadambat Nair

**Contact:** [bhklab.jamesbannon@gmail.com](mailto:bhklab.jamesbannon@gmail.com)

**Description:** Pipeline to build the Harmonized Drug Dataset Version 3 (HDD_v3) as a MultiAssayExperiment that harmonizes drug measurements, annotations, and fingerprints across multiple sources.

--------------------------------------

[![pixi-badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json&style=flat-square)](https://github.com/prefix-dev/pixi)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json&style=flat-square)](https://github.com/astral-sh/ruff)
[![Built with Material for MkDocs](https://img.shields.io/badge/mkdocs--material-gray?logo=materialformkdocs&style=flat-square)](https://github.com/squidfunk/mkdocs-material)

![GitHub last commit](https://img.shields.io/github/last-commit/bhklab/hdd-data-pipeline?style=flat-square)
![GitHub issues](https://img.shields.io/github/issues/bhklab/hdd-data-pipeline?style=flat-square)
![GitHub pull requests](https://img.shields.io/github/issues-pr/bhklab/hdd-data-pipeline?style=flat-square)
![GitHub contributors](https://img.shields.io/github/contributors/bhklab/hdd-data-pipeline?style=flat-square)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/bhklab/hdd-data-pipeline?style=flat-square)

## What the pipeline produces

- `data/results/HDD_v3.RDS`: the Harmonized Drug Dataset Version 3 as a `MultiAssayExperiment`.
- `data/results/HDD_v3_tables/`: MAE-derived TSV exports for metadata and dense assays, plus Matrix Market exports for sparse Morgan fingerprints.
- `data/results/HDD_v3_tables.tar.gz`: max-compression archive of the flat-table exports.
- `data/procdata/colData.tsv`: compound metadata assembled from AnnotationDB, curated JUMP-CP, OASIS, GEOM, and LINCS MAE inputs, and the CTRPv2 and NCI60 PharmacoSets, including ATC codes and LTKB, LiverTox, DIRIL, and DICT toxicity labels when available. ChEMBL drug indications are preserved separately as long-form metadata.
- `data/procdata/experiments/`: assay matrices for bioassays, DeepChem tasks, and sparse fingerprint features.
- `qc/hdd_quality_control.html`: quality control report rendered from `qc/hdd_quality_control.Rmd`.

## Quickstart

### Prerequisites

Pixi is required to run this project. If you have not installed it yet, follow the instructions at https://pixi.sh/latest/.

### Install dependencies

```bash
pixi install
```

### Run the pipeline

```bash
pixi run pipeline
```

Increase `-c` for more cores. The default Snakemake target builds and validates the RDS, table export, and archive. ClinTox and ToxCast are deprecated in v3 and are excluded from every release artifact.

### Generate the QC report

```bash
pixi run qc
```

## Configuration

Data sources, versions, and filtering rules are controlled in `config/pipeline.yaml`. Update this file to:

- Pin different dataset versions or URLs.
- Update curated sub-dataset MAE URLs.
- Change Morgan fingerprint radii and dimensions.
- Configure AnnotationDB batching, indication retrieval, and optional authenticated ATC enrichment.

## Repository layout

- `config/`: pipeline configuration.
- `workflow/`: Snakemake rules and scripts for data preparation and assembly.
- `data/rawdata/`: raw downloads (not tracked in Git).
- `data/procdata/`: processed intermediate datasets (not tracked in Git).
- `data/results/`: final HDD_v3 output (not tracked in Git).
- `data/results/HDD_v3_tables/`: TSV and Matrix Market exports with a SHA-256 file manifest.
- `qc/`: QC source notebook and rendered report.
- `docs/`: project documentation (this file, usage notes, data sources, dev notes).

## Additional documentation

- `docs/usage.md`: how to configure and run the pipeline.
- `docs/data_sources.md`: data source registry for HDD_v3 inputs.
- `docs/devnotes.md`: engineering notes and decisions.
