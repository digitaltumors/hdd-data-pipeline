from pathlib import Path

from damply import dirs

configfile: "config/pipeline.yaml"
SCRIPT_DIR = (
	Path(workflow.current_basedir.get_path_or_uri(secret_free=True))
	/ "workflow"
	/ "scripts"
)

include: "workflow/rules/fetchDeepChem.smk"
include: "workflow/rules/sub_dataset.smk"
include: "workflow/rules/annotationdb.smk"
include: "workflow/rules/deepchem_experiments.smk"
include: "workflow/rules/fingerprints.smk"
include: "workflow/rules/construct_mae.smk"

rule all:
	input:
		dirs.RESULTS / "HDD_v2.1.RDS",
		dirs.RESULTS / "HDD_v2.1_csv"
