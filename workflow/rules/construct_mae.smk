rule construct_MAE:
    input:
        colData=rules.process_AnnotationDB.output.colData,
        bioassays=rules.process_AnnotationDB.output.bioassays,
        tox21=rules.make_deepchem_experiments.output.tox21,
        sider=rules.make_deepchem_experiments.output.sider,
        indications=rules.process_AnnotationDB.output.indications,
        fetch_manifest=rules.fetch_AnnotationDB_raw.output.manifest,
        configfile="config/pipeline.yaml",
        fingerprints=rules.make_fingerprints.output.fingerprints,
        fingerprint_columns=rules.make_fingerprints.output.fingerprint_columns,
    output:
        mae=MAE_RDS,
    params:
        dataset_id=DATASET_ID,
        dataset_version=DATASET_VERSION,
    script:
        str(SCRIPT_DIR / "construct_MAE.R")


rule export_MAE_tables:
    input:
        mae=rules.construct_MAE.output.mae,
    output:
        outdir=directory(TABLE_DIR),
    script:
        str(SCRIPT_DIR / "export_mae_tables.R")


rule archive_MAE_tables:
    input:
        outdir=rules.export_MAE_tables.output.outdir,
    output:
        archive=TABLE_ARCHIVE,
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.archive}")"
        tar -C "{RESULTS}" -cf - "$(basename "{input.outdir}")" | gzip -9 >"{output.archive}"
        """


rule validate_release:
    input:
        mae=rules.construct_MAE.output.mae,
    output:
        validation=RELEASE_VALIDATION,
    params:
        dataset_id=DATASET_ID,
        dataset_version=DATASET_VERSION,
        expected_indication_count=dataset_config.get("expected_indication_count"),
        deprecated_experiments=dataset_config.get("deprecated_experiments", []),
        deepchem_expected=config["deep_chem"]["expected_assays"],
    script:
        str(SCRIPT_DIR / "validate_release.R")
