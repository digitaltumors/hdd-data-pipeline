deep_chem_urls = config["deep_chem"]["urls"]
deepchem_subdir = config["deep_chem"]["subdir"]
deepchem_datasets = tuple(deep_chem_urls.keys())
deepchem_dataset_pattern = "|".join(deepchem_datasets)


rule download_DeepChem_dataset:
    output:
        data=RAWDATA / deepchem_subdir / "{dataset}.csv",
    wildcard_constraints:
        dataset=deepchem_dataset_pattern,
    threads: 1
    params:
        url=lambda wildcards: deep_chem_urls[wildcards.dataset],
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.data}")"
        if [[ "{params.url}" == *.gz ]]; then
            curl -L --fail --silent --show-error "{params.url}" | gzip -dc >"{output.data}"
        else
            curl -L --fail --silent --show-error "{params.url}" -o "{output.data}"
        fi
        """
