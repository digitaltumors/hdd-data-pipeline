from damply import dirs


def membership_output(source_name):
	info = config[source_name]
	return (
		dirs.RAWDATA
		/ info["subdir"]
		/ info["version"]
		/ f"{source_name}_hdd_membership.tsv.gz"
	)


rule download_OASIS_membership:
	params:
		url = config["oasis"]["url"],
		sha256 = config["oasis"].get("sha256", "")

	output:
		data = membership_output("oasis")

	threads: 1

	shell:
		"""
		set -euo pipefail
		mkdir -p "$(dirname "{output.data}")"
		tmp="{output.data}.tmp"
		downloaded="${{tmp}}.download"
		curl -L --fail --silent --show-error "{params.url}" -o "$downloaded"
		if gzip -t "$downloaded" >/dev/null 2>&1; then
			mv "$downloaded" "$tmp"
		else
			gzip -c "$downloaded" > "$tmp"
			rm -f "$downloaded"
		fi
		if [[ -n "{params.sha256}" ]]; then
			printf '%s  %s\n' "{params.sha256}" "$tmp" | shasum -a 256 -c -
		fi
		mv "$tmp" "{output.data}"
		"""


rule download_GEOM_membership:
	params:
		url = config["geom"]["url"],
		sha256 = config["geom"].get("sha256", "")

	output:
		data = membership_output("geom")

	threads: 1

	shell:
		"""
		set -euo pipefail
		mkdir -p "$(dirname "{output.data}")"
		tmp="{output.data}.tmp"
		downloaded="${{tmp}}.download"
		curl -L --fail --silent --show-error "{params.url}" -o "$downloaded"
		if gzip -t "$downloaded" >/dev/null 2>&1; then
			mv "$downloaded" "$tmp"
		else
			gzip -c "$downloaded" > "$tmp"
			rm -f "$downloaded"
		fi
		if [[ -n "{params.sha256}" ]]; then
			printf '%s  %s\n' "{params.sha256}" "$tmp" | shasum -a 256 -c -
		fi
		mv "$tmp" "{output.data}"
		"""
