from damply import dirs

deepchem_subdir = config["deep_chem"]["subdir"]
deepchem_experiments = ["toxcast", "tox21", "sider", "clintox"]


rule make_deepchem_experiments:
	params:
		deepchem_subdir = deepchem_subdir

	input:
		colData = rules.process_AnnotationDB.output.colData,
		smiles_to_cid = rules.download_DeepChem_smiles_cid_mapping.output.data,
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
		str(SCRIPT_DIR / "make_deepchem_experiments.py")
