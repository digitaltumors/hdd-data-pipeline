from damply import dirs

jumpcp_subdir = config['jump_cp']['subdir']
jumpcp_version = config['jump_cp']['version']


rule download_JUMPCP:
	params:
		cpd_url = config['jump_cp']['url']

	output:
		data = dirs.RAWDATA / jumpcp_subdir / jumpcp_version / 'JUMP_CP_compounds.csv'

	threads: 1

	shell:
		"""
		set -euo pipefail
		mkdir -p "$(dirname "{output.data}")"
		curl -L --fail --silent --show-error "{params.cpd_url}" | gzip -dc > "{output.data}"
		"""
