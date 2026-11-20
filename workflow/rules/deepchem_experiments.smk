deepchem_subdir = config["deep_chem"]["subdir"]


rule make_deepchem_experiments:
    input:
        colData=rules.process_AnnotationDB.output.colData,
        sider=RAWDATA / deepchem_subdir / "sider.csv",
    output:
        sider=PROCDATA / "experiments" / "sider.tsv",
        match_audit=PROCDATA / "metadata" / "deepchem_match_audit.tsv",
        candidate_audit=PROCDATA / "metadata" / "deepchem_match_candidates.tsv",
        collision_audit=PROCDATA / "metadata" / "deepchem_match_collisions.tsv",
        match_summary=PROCDATA / "metadata" / "deepchem_match_summary.tsv",
    params:
        expected_assays=config["deep_chem"]["expected_assays"],
    script:
        str(SCRIPT_DIR / "make_deepchem_experiments.py")
