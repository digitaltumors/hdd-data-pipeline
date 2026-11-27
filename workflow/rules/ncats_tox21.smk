ncats = config["ncats_tox21"]


rule download_NCATS_Tox21:
    output:
        archive=RAWDATA / ncats["subdir"] / "tox21_10k_data_all.zip",
    params:
        url=ncats["url"],
    shell:
        """
        mkdir -p "$(dirname {output.archive:q})"
        curl --location --fail --silent --show-error --retry 3 {params.url:q} -o {output.archive:q}
        """


rule make_ncats_tox21:
    input:
        archive=rules.download_NCATS_Tox21.output.archive,
        coldata=rules.process_AnnotationDB.output.colData,
        implementation=SCRIPT_DIR / "make_ncats_tox21.py",
    output:
        assay=PROCDATA / "experiments" / "tox21.tsv",
    shell:
        """
        python {input.implementation:q} --archive {input.archive:q} --coldata {input.coldata:q} --assay {output.assay:q}
        """
