from damply import dirs


rule construct_MAE:
	input:
		colData = rules.process_AnnotationDB.output.colData,
		bioassays = rules.process_AnnotationDB.output.bioassays,
		bindingdb = rules.make_bindingdb_experiments.output.binding_db,
		toxcast = rules.make_deepchem_experiments.output.toxcast,
		tox21 = rules.make_deepchem_experiments.output.tox21,
		sider = rules.make_deepchem_experiments.output.sider,
		clintox = rules.make_deepchem_experiments.output.clintox,
		fingerprints = rules.make_fingerprints.output.fingerprints

	output:
		mae = dirs.RESULTS / "HDD_v1.1.RDS"

	script:
		str(SCRIPT_DIR / "construct_MAE.R")


rule export_MAE_csvs:
	input:
		mae = rules.construct_MAE.output.mae

	output:
		outdir = directory(dirs.RESULTS / "HDD_v1.1_csv")

	script:
		str(SCRIPT_DIR / "export_mae_csvs.R")
