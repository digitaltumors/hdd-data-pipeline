deepchem_subdir = config["deep_chem"]["subdir"]
deepchem_experiments = ["tox21", "sider"]


rule make_deepchem_experiments:
    input:
        colData=rules.process_AnnotationDB.output.colData,
        tox21=RAWDATA / deepchem_subdir / "tox21.csv",
        sider=RAWDATA / deepchem_subdir / "sider.csv",
    output:
        tox21=PROCDATA / "experiments" / "tox21.tsv",
        sider=PROCDATA / "experiments" / "sider.tsv",
        match_audit=PROCDATA / "metadata" / "deepchem_match_audit.tsv",
        candidate_audit=PROCDATA / "metadata" / "deepchem_match_candidates.tsv",
        collision_audit=PROCDATA / "metadata" / "deepchem_match_collisions.tsv",
        match_summary=PROCDATA / "metadata" / "deepchem_match_summary.tsv",
    params:
        deepchem_subdir=deepchem_subdir,
        expected_assays=config["deep_chem"]["expected_assays"],
    script:
        str(SCRIPT_DIR / "make_deepchem_experiments.py")
