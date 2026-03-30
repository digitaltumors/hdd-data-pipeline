from damply import dirs


rule make_bindingdb_experiments:
	input:
		colData = rules.process_AnnotationDB.output.colData,
		bdb_data = rules.process_BindingDB.output.cleaned_data

	output:
		binding_db = dirs.PROCDATA / "experiments" / "binding_db.csv"

	script:
		str(SCRIPT_DIR / "make_bindingdbd_experiments.py")
