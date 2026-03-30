from pathlib import Path

from damply import dirs

configfile: "config/pipeline.yaml"
SCRIPT_DIR = Path(workflow.source_path("workflow/scripts/fetch_annotationdb.py")).parent

include: "workflow/rules/fetchLINCS.smk"
include: "workflow/rules/fetchJUMPCP.smk"
include: "workflow/rules/fetchDeepChem.smk"
include: "workflow/rules/annotationdb.smk"
include: "workflow/rules/deepchem_experiments.smk"
include: "workflow/rules/fingerprints.smk"
include: "workflow/rules/construct_mae.smk"

rule all:
	input:
		dirs.RESULTS / "HDD_v1.RDS",
		dirs.RESULTS / "HDD_v1_csv"
