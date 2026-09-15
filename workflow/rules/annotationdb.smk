import json

annotationdb_config = config["annotationdb"]
annotationdb_fetch = annotationdb_config.get("fetch", {})
annotationdb_base_url = annotationdb_config["base_url"].rstrip("/")
annotationdb_batch_size = annotationdb_fetch.get("batch_size", 50)
annotationdb_workers = annotationdb_fetch.get("workers", 4)
annotationdb_db_url = annotationdb_base_url + annotationdb_config.get(
    "compound_all_path", "/compound/all"
)
annotationdb_details_url = annotationdb_base_url + annotationdb_config.get(
    "compound_many_path", "/compound/many"
)
annotationdb_query_delay_seconds = annotationdb_fetch.get("query_delay_seconds", 0.1)
annotationdb_golden_bioassay = annotationdb_fetch.get("golden_bioassay", True)
annotationdb_indication = annotationdb_fetch.get("indication", True)
annotationdb_api_key_file = annotationdb_fetch.get("api_key_file", ".env")
annotationdb_api_key_name = annotationdb_fetch.get(
    "api_key_name", "ANNOTATIONDB_API_KEY"
)


rule fetch_AnnotationDB_raw:
    output:
        raw=PROCDATA / "ANNOTATION_DB" / "compound_details.jsonl",
        manifest=PROCDATA / "ANNOTATION_DB" / "fetch_manifest.tsv",
    threads: annotationdb_workers
    params:
        db_url=annotationdb_db_url,
        details_url=annotationdb_details_url,
        batch_size=annotationdb_batch_size,
        query_delay_seconds=annotationdb_query_delay_seconds,
        golden_bioassay=annotationdb_golden_bioassay,
        indication=annotationdb_indication,
        api_key_file=annotationdb_api_key_file,
        api_key_name=annotationdb_api_key_name,
    script:
        str(SCRIPT_DIR / "fetch_annotationdb.py")


rule process_AnnotationDB:
    input:
        raw_data=rules.fetch_AnnotationDB_raw.output.raw,
        sub_dataset_metadata=expand(
            PROCDATA / "sub_dataset" / "{dataset}_drug_metadata.tsv",
            dataset=sub_dataset_names,
        ),
        bbbp_file=RAWDATA
        / config["deep_chem"]["subdir"]
        / "blood_brain_barrier.csv",
    output:
        colData=PROCDATA / "colData.tsv",
        bioassays=PROCDATA / "experiments" / "bioassays.tsv",
        indications=PROCDATA / "metadata" / "drug_indications.tsv",
        parity=directory(PROCDATA / "sub_dataset" / "parity"),
    threads: 1
    params:
        db_url=annotationdb_db_url,
        expected_indication_count=dataset_config.get("expected_indication_count"),
        sub_dataset_names=list(sub_dataset_names),
        sub_dataset_specs=json.dumps(config["sub_dataset"], sort_keys=True),
    script:
        str(SCRIPT_DIR / "process_annotationdb.py")
