from damply import dirs

lincs_version = config['lincs']['version']
lincs_subdir = config['lincs']['subdir']


rule download_LINCS:
	params:
		download_url = config['lincs']['url']

	output:
		lincs_raw = dirs.RAWDATA / lincs_subdir / lincs_version / 'compounds_raw.tsv'

	threads: 1

	shell:
		"""
		set -euo pipefail
		mkdir -p "$(dirname "{output.lincs_raw}")"
		curl -L --fail --silent --show-error "{params.download_url}" -o "{output.lincs_raw}"
		"""
