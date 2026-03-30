from damply import dirs

fingerprint_radii = config["colData"]["fingerprints"]["radius_list"]
fingerprint_dims = config["colData"]["fingerprints"]["dim_list"]


rule make_fingerprints:
	params:
		radius_list = fingerprint_radii,
		dim_list = fingerprint_dims

	input:
		rules.fetch_from_AnnotationDB.output.colData

	output:
		fingerprints = expand(
			dirs.PROCDATA / "experiments" / "fingerprints" / "Morgan.{rad}.{dim}.mtx",
			rad=fingerprint_radii,
			dim=fingerprint_dims,
		)

	script:
		str(SCRIPT_DIR / "make_fingerprints.py")
