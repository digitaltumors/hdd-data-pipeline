# Harmonized Drug Dataset v3 Pipeline

**Authors:** [James Bannon](https://github.com/jbannon), Michael Tran, Matthew Boccalon, Sisira Kadambat Nair

**Contact:** [bhklab.jamesbannon@gmail.com](mailto:bhklab.jamesbannon@gmail.com)

[![pixi-badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json&style=flat-square)](https://github.com/prefix-dev/pixi)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json&style=flat-square)](https://github.com/astral-sh/ruff)
[![Built with Material for MkDocs](https://img.shields.io/badge/mkdocs--material-gray?logo=materialformkdocs&style=flat-square)](https://github.com/squidfunk/mkdocs-material)

![GitHub last commit](https://img.shields.io/github/last-commit/bhklab/hdd-data-pipeline?style=flat-square)
![GitHub issues](https://img.shields.io/github/issues/bhklab/hdd-data-pipeline?style=flat-square)
![GitHub pull requests](https://img.shields.io/github/issues-pr/bhklab/hdd-data-pipeline?style=flat-square)
![GitHub contributors](https://img.shields.io/github/contributors/bhklab/hdd-data-pipeline?style=flat-square)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/bhklab/hdd-data-pipeline?style=flat-square)

## Overview

This repository regenerates HDD v3 as a compound-keyed
`MultiAssayExperiment`. It harmonizes AnnotationDB compounds with pinned JUMP,
OASIS, GEOM, LINCS, CTRPv2, and NCI60 objects, adds ChEMBL indication records,
constructs selected assay matrices, and computes sparse Morgan fingerprints.

## Quick Start

Pixi is required to run this project. If it is not installed, follow the
instructions at <https://pixi.sh/latest/>.

```bash
pixi install
pixi run pipeline
pixi run qc
```

The full regeneration downloads about 7.4 GB of pinned source objects plus the
live AnnotationDB snapshot. Expected end-to-end runtime is roughly 1-4 hours,
depending on network throughput and RDKit fingerprint generation.

Increase Snakemake's core count when appropriate by passing arguments after the
task name, for example `pixi run pipeline -- -c 4`.

## Pipeline

The Snakemake workflow fetches and validates sources, extracts sub-dataset drug
metadata, builds stable `HDD.Compound.ID` values, preserves source membership
and IDs, writes ChEMBL indications as a long table, constructs SIDER,
Bioassays, and Tox21 experiments, generates four sparse Morgan assays, builds
the MAE, validates the release contract, and exports TSV/Matrix Market files.

ClinTox and ToxCast are deprecated in v3 and are absent from the download DAG,
processed intermediates, final object, exports, QC, and release allowlist.

## Inputs

`config/pipeline.yaml` is the source of truth for paths, release identity,
source URLs, AnnotationDB request controls, expected indication count, assay
deprecations, and fingerprint dimensions. Optional ATC authorization can be
provided through the ignored `.env` file; all other annotations remain
regenerative without it.

The configured inputs include curated JUMP-CP, OASIS, GEOM, and LINCS MAEs;
CTRPv2 and NCI60 PharmacoSets; DeepChem BBBP, Tox21, and SIDER tables; and the
AnnotationDB compound index and detail endpoints. Source versions, URLs, and
access notes are listed in `docs/data_sources.md`.

## Outputs

- `data/results/HDD_v3.RDS`
- `data/results/HDD_v3_tables/`
- `data/results/HDD_v3_tables.tar.gz`
- `qc/hdd_quality_control.html`

Dense assays are exported as TSV. Morgan fingerprints remain sparse Matrix
Market files with explicit row and compound maps. Every export has a SHA-256
entry in `file_manifest.tsv`; the embedded release contract is also available
as `metadata/pipeline_config.tsv`.

## QC and Docs

QC source lives at `qc/hdd_quality_control.Rmd`; `pixi run qc` renders the HTML
report beside it.

## Development

Use `pixi run pipeline -- --dry-run` to inspect the DAG. Data-source and schema
decisions are documented under `docs/`.

Repository layout:

- `config/`: release identity, source URLs, expected counts, and parameters.
- `workflow/rules/`: Snakemake stages.
- `workflow/scripts/`: processing, assembly, validation, and export code.
- `qc/`: QC source notebook and rendered report.
- `data/rawdata/`: downloaded source data, ignored by Git.
- `data/procdata/`: generated intermediate data and audit tables, ignored by Git.
- `data/results/`: final release object and exports, ignored by Git.
- `docs/`: usage, provenance, and development decisions.
