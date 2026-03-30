from damply import dirs

annotationdb_fetch = config["colData"].get("fetch", {})
annotationdb_batch_size = annotationdb_fetch.get("batch_size", 50)
annotationdb_workers = annotationdb_fetch.get("workers", 4)
annotationdb_details_url = annotationdb_fetch.get("details_url")
annotationdb_query_delay_seconds = annotationdb_fetch.get("query_delay_seconds", 0.1)
annotationdb_golden_bioassay = annotationdb_fetch.get("golden_bioassay", True)


rule fetch_AnnotationDB_raw:
	params:
		db_url = config["colData"]["db_url"],
		details_url = annotationdb_details_url,
		batch_size = annotationdb_batch_size,
		query_delay_seconds = annotationdb_query_delay_seconds,
		golden_bioassay = annotationdb_golden_bioassay

	output:
		raw = dirs.PROCDATA / "ANNOTATION_DB" / "compound_details.jsonl"

	threads: annotationdb_workers

	script:
		str(SCRIPT_DIR / "fetch_annotationdb.py")


rule fetch_from_AnnotationDB:
	params:
		db_url = config["colData"]["db_url"]

	input:
		raw_data = rules.fetch_AnnotationDB_raw.output.raw,
		lincs_file = rules.download_LINCS.output.lincs_raw,
		jump_file = rules.download_JUMPCP.output.data,
		bbbp_file = dirs.RAWDATA / config["deep_chem"]["subdir"] / "blood_brain_barrier.csv"

	output:
		colData = dirs.PROCDATA / "colData.csv",
		bioassays = dirs.PROCDATA / "experiments" / "bioassays.csv"

	threads: 1

	script:
		str(SCRIPT_DIR / "process_annotationdb.py")
