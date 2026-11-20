suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(S4Vectors)
})

stopf <- function(message, ...) {
  stop(sprintf(message, ...), call. = FALSE)
}

if (!exists("snakemake")) {
  stopf("validate_release.R must be run through Snakemake")
}

mae <- readRDS(snakemake@input[["mae"]])
output_path <- snakemake@output[["validation"]]
params <- snakemake@params
object_metadata <- S4Vectors::metadata(mae)
experiment_names <- names(MultiAssayExperiment::experiments(mae))

expected_experiments <- c(
  "SIDER",
  "Bioassays",
  "Tox21",
  "fingerprint.Morgan.2.1024",
  "fingerprint.Morgan.3.1024",
  "fingerprint.Morgan.2.2048",
  "fingerprint.Morgan.3.2048"
)
if (!setequal(experiment_names, expected_experiments)) {
  stopf(
    "Experiment set mismatch: %s",
    paste(experiment_names, collapse = ", ")
  )
}
deprecated <- as.character(params[["deprecated_experiments"]])
if (length(intersect(experiment_names, deprecated)) > 0L) {
  stopf("Deprecated experiments remain in HDD v3")
}

sider_expected <- params[["deepchem_expected"]][["sider"]]
sider_dimensions <- dim(MultiAssayExperiment::experiments(mae)[["SIDER"]])
if (!identical(sider_dimensions, as.integer(c(
  sider_expected[["features"]], sider_expected[["compounds"]]
)))) {
  stopf(
    "SIDER dimensions mismatch: observed %d x %d, expected %d x %d",
    sider_dimensions[[1]], sider_dimensions[[2]],
    sider_expected[["features"]], sider_expected[["compounds"]]
  )
}

tox21 <- as.matrix(SummarizedExperiment::assay(
  MultiAssayExperiment::experiments(mae)[["Tox21"]]
))
tox21_source <- utils::read.delim(
  snakemake@input[["tox21_assay"]], row.names = 1L, check.names = FALSE,
  na.strings = c("NA", ""), colClasses = c("character", rep("integer", ncol(tox21)))
)
if (
  nrow(tox21) != as.integer(params[["tox21_features"]]) ||
    ncol(tox21) == 0L ||
    any(rowSums(!is.na(tox21)) == 0L) ||
    !identical(dimnames(tox21), dimnames(as.matrix(tox21_source))) ||
    !isTRUE(all.equal(unname(tox21), unname(as.matrix(tox21_source)))) ||
    any(!is.na(tox21) & !tox21 %in% c(0, 1))
) {
  stopf("NCATS Tox21 assay dimensions or calls do not match the generated data")
}
if (!identical(object_metadata$Pipeline$ID, params[["dataset_id"]])) {
  stopf("Dataset ID mismatch")
}
if (!identical(object_metadata$Pipeline$Version, params[["dataset_version"]])) {
  stopf("Dataset version mismatch")
}

coldata <- as.data.frame(MultiAssayExperiment::colData(mae))
if (
  any(is.na(coldata$HDD.Compound.ID)) ||
    anyDuplicated(coldata$HDD.Compound.ID) > 0L ||
    !identical(row.names(coldata), as.character(coldata$HDD.Compound.ID))
) {
  stopf("HDD.Compound.ID must be unique, non-missing, and equal row names")
}

indications <- as.data.frame(object_metadata$Drug.Indications)
expected_indications <- as.integer(params[["expected_indication_count"]])
if (nrow(indications) != expected_indications) {
  stopf(
    "Drug indication count mismatch: observed %d, expected %d",
    nrow(indications),
    expected_indications
  )
}
if (
  any(is.na(indications$HDD.Compound.ID)) ||
    length(setdiff(indications$HDD.Compound.ID, coldata$HDD.Compound.ID)) > 0L
) {
  stopf("Drug.Indications contains invalid HDD.Compound.ID foreign keys")
}

sample_map <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
if (length(setdiff(sample_map$primary, coldata$HDD.Compound.ID)) > 0L) {
  stopf("sampleMap contains IDs absent from colData")
}
stopifnot(validObject(mae))

validation <- data.frame(
  Check = c(
    "valid_object",
    "dataset_identity",
    "experiment_allowlist",
    "deprecated_experiments_absent",
    "sider_dimensions",
    "ncats_tox21_assay",
    "compound_ids",
    "drug_indications",
    "sample_map"
  ),
  Status = "PASS",
  stringsAsFactors = FALSE
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  validation,
  output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
