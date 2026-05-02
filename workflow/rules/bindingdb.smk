from damply import dirs


bindingdb_config = config["binding_db"]
bindingdb_subset = bindingdb_config["subset"]
bindingdb_version = bindingdb_config["version"]
bindingdb_stem = f"BindingDB_{bindingdb_subset}_{bindingdb_version}"


rule download_BindingDB:
	params:
		url = (
			f"{bindingdb_config['base_url']}_"
			f"{bindingdb_subset}_{bindingdb_version}_tsv.zip"
		)

	output:
		zipped_data = dirs.RAWDATA / "BINDING_DB" / f"{bindingdb_stem}_tsv.zip"

	shell:
		"""
		mkdir -p $(dirname {output.zipped_data})
		curl -L --fail --silent --show-error {params.url} -o {output.zipped_data}
		"""


rule process_BindingDB:
	params:
		processing = bindingdb_config["processing"]

	input:
		raw_zip = rules.download_BindingDB.output.zipped_data

	output:
		cleaned_data = dirs.PROCDATA / "BINDING_DB" / f"{bindingdb_stem}_cleaned.csv"

	script:
		str(SCRIPT_DIR / "process_bindingdb.py")


rule make_bindingdb_experiments:
	input:
		colData = rules.process_AnnotationDB.output.colData,
		bindingdb = rules.process_BindingDB.output.cleaned_data

	output:
		binding_db = dirs.PROCDATA / "experiments" / "binding_db.csv"

	script:
		str(SCRIPT_DIR / "make_bindingdb_experiments.py")
