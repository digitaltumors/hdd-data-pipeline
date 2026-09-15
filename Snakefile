configfile: "config/pipeline.yaml"


include: "workflow/Snakefile"


rule all:
    default_target: True
    input:
        ALL_TARGETS,
