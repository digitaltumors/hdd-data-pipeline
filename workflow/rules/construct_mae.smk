from damply import dirs


rule construct_MAE:
	input:
		colData = rules.fetch_from_AnnotationDB.output.colData,
		bioassays = rules.fetch_from_AnnotationDB.output.bioassays,
		toxcast = rules.make_deepchem_experiments.output.toxcast,
		tox21 = rules.make_deepchem_experiments.output.tox21,
		sider = rules.make_deepchem_experiments.output.sider,
		clintox = rules.make_deepchem_experiments.output.clintox,
		fingerprints = rules.make_fingerprints.output.fingerprints

	output:
		mae = dirs.RESULTS / "HDD_v1.RDS"

	script:
		"workflow/scripts/construct_MAE.R"


rule export_MAE_csvs:
	input:
		mae = rules.construct_MAE.output.mae

	output:
		outdir = directory(dirs.RESULTS / "HDD_v1_csv")

	script:
		"workflow/scripts/export_mae_csvs.R"
