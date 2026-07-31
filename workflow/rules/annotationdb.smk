from damply import dirs
import json

annotationdb_fetch = config["colData"].get("fetch", {})
annotationdb_batch_size = annotationdb_fetch.get("batch_size", 50)
annotationdb_workers = annotationdb_fetch.get("workers", 4)
annotationdb_details_url = annotationdb_fetch.get("details_url")
annotationdb_query_delay_seconds = annotationdb_fetch.get("query_delay_seconds", 0.1)
annotationdb_golden_bioassay = annotationdb_fetch.get("golden_bioassay", True)
annotationdb_env_file = annotationdb_fetch.get("env_file", ".env")
annotationdb_api_key_env_var = annotationdb_fetch.get(
    "api_key_env_var", "ANNOTATIONDB_API_KEY"
)


rule fetch_AnnotationDB_raw:
    output:
        raw=dirs.PROCDATA / "ANNOTATION_DB" / "compound_details.jsonl",
    threads: annotationdb_workers
    params:
        db_url=config["colData"]["db_url"],
        details_url=annotationdb_details_url,
        batch_size=annotationdb_batch_size,
        query_delay_seconds=annotationdb_query_delay_seconds,
        golden_bioassay=annotationdb_golden_bioassay,
        env_file=annotationdb_env_file,
        api_key_env_var=annotationdb_api_key_env_var,
    script:
        str(SCRIPT_DIR / "fetch_annotationdb.py")


rule process_AnnotationDB:
    input:
        raw_data=rules.fetch_AnnotationDB_raw.output.raw,
        sub_dataset_metadata=expand(
            dirs.PROCDATA / "sub_dataset" / "{dataset}_drug_metadata.tsv",
            dataset=sub_dataset_names,
        ),
        bbbp_file=dirs.RAWDATA
        / config["deep_chem"]["subdir"]
        / "blood_brain_barrier.csv",
    output:
        colData=dirs.PROCDATA / "colData.csv",
        bioassays=dirs.PROCDATA / "experiments" / "bioassays.csv",
        parity=directory(dirs.PROCDATA / "sub_dataset" / "parity"),
    threads: 1
    params:
        db_url=config["colData"]["db_url"],
        sub_dataset_names=list(sub_dataset_names),
        sub_dataset_specs=json.dumps(config["sub_dataset"], sort_keys=True),
    script:
        str(SCRIPT_DIR / "process_annotationdb.py")
