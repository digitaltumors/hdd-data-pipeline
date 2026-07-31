from damply import dirs

sub_dataset_config = config["sub_dataset"]
sub_dataset_names = tuple(sub_dataset_config.keys())
sub_dataset_pattern = "|".join(sub_dataset_names)


rule download_sub_dataset_object:
    output:
        rds=dirs.RAWDATA / "sub_dataset" / "{dataset}" / "{dataset}.rds",
    wildcard_constraints:
        dataset=sub_dataset_pattern,
    threads: 1
    params:
        url=lambda wildcards: sub_dataset_config[wildcards.dataset]["url"],
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.rds}")"
        curl -L --fail --silent --show-error "{params.url}" -o "{output.rds}"
        """


rule extract_sub_dataset_drug_metadata:
    input:
        rds=rules.download_sub_dataset_object.output.rds,
    output:
        metadata=dirs.PROCDATA / "sub_dataset" / "{dataset}_drug_metadata.tsv",
        summary=dirs.PROCDATA / "sub_dataset" / "{dataset}_drug_metadata_summary.tsv",
    wildcard_constraints:
        dataset=sub_dataset_pattern,
    params:
        dataset=lambda wildcards: wildcards.dataset,
        object_type=lambda wildcards: sub_dataset_config[wildcards.dataset][
            "object_type"
        ],
        source_key_column=lambda wildcards: sub_dataset_config[wildcards.dataset][
            "source_key_column"
        ],
        source_key_source_columns=lambda wildcards: sub_dataset_config[
            wildcards.dataset
        ].get("source_key_source_columns", []),
    script:
        str(SCRIPT_DIR / "extract_sub_dataset_drug_metadata.R")
