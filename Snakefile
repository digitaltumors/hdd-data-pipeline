from pathlib import Path

from damply import dirs

configfile: "config/pipeline.yaml"
SCRIPT_DIR = (
	Path(workflow.current_basedir.get_path_or_uri(secret_free=True))
	/ "workflow"
	/ "scripts"
)

include: "workflow/rules/fetchLINCS.smk"
include: "workflow/rules/fetchJUMPCP.smk"
include: "workflow/rules/fetchMembership.smk"
include: "workflow/rules/fetchDeepChem.smk"
include: "workflow/rules/annotationdb.smk"
include: "workflow/rules/deepchem_experiments.smk"
include: "workflow/rules/fingerprints.smk"
include: "workflow/rules/construct_mae.smk"

rule all:
	input:
		dirs.RESULTS / "HDD_v2.RDS",
		dirs.RESULTS / "HDD_v2_csv"
