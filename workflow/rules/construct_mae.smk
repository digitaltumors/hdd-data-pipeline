from damply import dirs


rule construct_MAE:
	input:
		colData = rules.fetch_from_AnnotationDB.output.colData,
		bioassays = rules.fetch_from_AnnotationDB.output.bioassays,
		deepchem = rules.make_deepchem_experiments.output.experiments,
		fingerprints = rules.make_fingerprints.output.fingerprints

	output:
		mae = dirs.RESULTS / "HDD_v1.RDS"

	shell:
		"""
		mkdir -p {dirs.RESULTS}
		Rscript ./workflow/scripts/construct_MAE.R
		"""


rule export_MAE_csvs:
	input:
		mae = rules.construct_MAE.output.mae

	output:
		outdir = directory(dirs.RESULTS / "HDD_v1_csv")

	shell:
		"""
		Rscript ./workflow/scripts/export_mae_csvs.R {input.mae} {output.outdir}
		"""
