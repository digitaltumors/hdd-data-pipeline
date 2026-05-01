from damply import dirs


rule construct_MAE:
	input:
		colData = rules.process_AnnotationDB.output.colData,
		bioassays = rules.process_AnnotationDB.output.bioassays,
		toxcast = rules.make_deepchem_experiments.output.toxcast,
		tox21 = rules.make_deepchem_experiments.output.tox21,
		sider = rules.make_deepchem_experiments.output.sider,
		clintox = rules.make_deepchem_experiments.output.clintox,
		fingerprints = rules.make_fingerprints.output.fingerprints

	output:
		mae = dirs.RESULTS / "HDD_v2.RDS"

	script:
		str(SCRIPT_DIR / "construct_MAE.R")


rule export_MAE_csvs:
	input:
		mae = rules.construct_MAE.output.mae

	output:
		outdir = directory(dirs.RESULTS / "HDD_v2_csv")

	script:
		str(SCRIPT_DIR / "export_mae_csvs.R")
