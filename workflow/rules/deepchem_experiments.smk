from damply import dirs

deepchem_subdir = config["deep_chem"]["subdir"]
deepchem_experiments = ["toxcast", "tox21", "sider", "clintox"]


rule make_deepchem_experiments:
	params:
		deepchem_subdir = deepchem_subdir

	input:
		colData = rules.fetch_from_AnnotationDB.output.colData,
		toxcast = dirs.RAWDATA / deepchem_subdir / "toxcast.csv",
		tox21 = dirs.RAWDATA / deepchem_subdir / "tox21.csv",
		sider = dirs.RAWDATA / deepchem_subdir / "sider.csv",
		clintox = dirs.RAWDATA / deepchem_subdir / "clintox.csv"

	output:
		toxcast = dirs.PROCDATA / "experiments" / "toxcast.csv",
		tox21 = dirs.PROCDATA / "experiments" / "tox21.csv",
		sider = dirs.PROCDATA / "experiments" / "sider.csv",
		clintox = dirs.PROCDATA / "experiments" / "clintox.csv"

	script:
		"workflow/scripts/make_deepchem_experiments.py"
