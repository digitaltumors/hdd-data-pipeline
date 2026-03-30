from damply import dirs

annotationdb_fetch = config["colData"].get("fetch", {})
annotationdb_batch_size = annotationdb_fetch.get("batch_size", 5)
annotationdb_workers = annotationdb_fetch.get("workers", 8)


rule fetch_AnnotationDB_raw:
	params:
		db_url = config["colData"]["db_url"],
		batch_size = annotationdb_batch_size

	output:
		raw = dirs.RAWDATA / "ANNOTATION_DB" / "compound_details.jsonl"

	threads: annotationdb_workers

	shell:
		"""
		mkdir -p {dirs.RAWDATA}/ANNOTATION_DB
		PYTHONPATH=workflow/scripts python ./workflow/scripts/fetch_annotationdb.py \
			-u {params.db_url} -o {output.raw} \
			--batch-size {params.batch_size} --workers {threads}
		"""


rule fetch_from_AnnotationDB:
	params:
		db_url = config["colData"]["db_url"]

	input:
		raw_data = rules.fetch_AnnotationDB_raw.output.raw,
		lincs_file = rules.download_LINCS.output.lincs_raw,
		jump_file = rules.download_JUMPCP.output.data,
		bbbp_file = dirs.PROCDATA / config["deep_chem"]["subdir"] / "blood_brain_barrier.csv"

	output:
		colData = dirs.PROCDATA / "colData.csv",
		bioassays = dirs.PROCDATA / "experiments" / "bioassays.csv"

	threads: 1

	shell:
		"""
		mkdir -p {dirs.PROCDATA} {dirs.PROCDATA}/experiments
		PYTHONPATH=workflow/scripts python ./workflow/scripts/process_annotationdb.py \
			-i {input.raw_data} -l {input.lincs_file} -j {input.jump_file} -b {input.bbbp_file}
		"""
