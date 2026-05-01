from damply import dirs


sub_dataset_config = config["sub_dataset"]
sub_dataset_names = tuple(sub_dataset_config.keys())
sub_dataset_pattern = "|".join(sub_dataset_names)


rule download_sub_dataset_mae:
	output:
		rds = dirs.RAWDATA / "sub_dataset" / "{dataset}" / "{dataset}_MultiAssayExperiment.rds"
	params:
		url = lambda wildcards: sub_dataset_config[wildcards.dataset]["url"]
	wildcard_constraints:
		dataset = sub_dataset_pattern
	threads: 1
	shell:
		"""
		set -euo pipefail
		mkdir -p "$(dirname "{output.rds}")"
		curl -L --fail --silent --show-error "{params.url}" -o "{output.rds}"
		"""


rule extract_sub_dataset_drug_metadata:
	input:
		rds = rules.download_sub_dataset_mae.output.rds
	output:
		metadata = dirs.PROCDATA / "sub_dataset" / "{dataset}_drug_metadata.tsv",
		summary = dirs.PROCDATA / "sub_dataset" / "{dataset}_drug_metadata_summary.tsv"
	params:
		dataset = lambda wildcards: wildcards.dataset
	wildcard_constraints:
		dataset = sub_dataset_pattern
	script:
		str(SCRIPT_DIR / "extract_sub_dataset_drug_metadata.R")
