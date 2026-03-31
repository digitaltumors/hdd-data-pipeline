from damply import dirs

deep_chem_urls = config['deep_chem']['urls']
deepchem_subdir = config['deep_chem']['subdir']
deepchem_smiles_to_cid_url = config['deep_chem']['pubchem_smiles_to_cid_url']
deepchem_datasets = tuple(deep_chem_urls.keys())
deepchem_dataset_pattern = '|'.join(deepchem_datasets)


rule download_DeepChem_dataset:
	output:
		data = dirs.RAWDATA / deepchem_subdir / '{dataset}.csv'

	params:
		url = lambda wildcards: deep_chem_urls[wildcards.dataset]

	wildcard_constraints:
		dataset = deepchem_dataset_pattern

	threads: 1

	shell:
		"""
		set -euo pipefail
		mkdir -p "$(dirname "{output.data}")"
		if [[ "{params.url}" == *.gz ]]; then
			curl -L --fail --silent --show-error "{params.url}" | gzip -dc > "{output.data}"
		else
			curl -L --fail --silent --show-error "{params.url}" -o "{output.data}"
		fi
		"""


rule download_DeepChem_smiles_cid_mapping:
	output:
		data = dirs.RAWDATA / deepchem_subdir / "deepchem_smiles_to_cid_pubchem.txt.gz"

	params:
		url = deepchem_smiles_to_cid_url

	threads: 1

	shell:
		"""
		set -euo pipefail
		mkdir -p "$(dirname "{output.data}")"
		curl -L --fail --silent --show-error "{params.url}" -o "{output.data}"
		"""
