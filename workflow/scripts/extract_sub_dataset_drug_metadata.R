suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(readr)
  library(tibble)
})

as_plain_data_frame <- function(x) {
  if (inherits(x, "DataFrame")) {
    x <- as.data.frame(x)
  } else if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }

  x
}

collapse_list_columns <- function(df) {
  for (column in names(df)) {
    if (!is.list(df[[column]])) {
      next
    }

    df[[column]] <- vapply(
      df[[column]],
      function(value) {
        if (length(value) == 0 || all(is.na(value))) {
          return(NA_character_)
        }
        paste(as.character(value), collapse = "|")
      },
      character(1)
    )
  }

  df
}

if (exists("snakemake")) {
  dataset <- snakemake@params[["dataset"]]
  rds_path <- snakemake@input[["rds"]]
  metadata_path <- snakemake@output[["metadata"]]
  summary_path <- snakemake@output[["summary"]]
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 4) {
    stop(
      "Usage: extract_sub_dataset_drug_metadata.R ",
      "<dataset> <input_rds> <metadata_tsv> <summary_tsv>",
      call. = FALSE
    )
  }
  dataset <- args[[1]]
  rds_path <- args[[2]]
  metadata_path <- args[[3]]
  summary_path <- args[[4]]
}

mae <- readRDS(rds_path)
mae_metadata <- metadata(mae)
if (!("Drug.Metadata" %in% names(mae_metadata))) {
  stop(
    "[extract_sub_dataset_drug_metadata] missing metadata(mae)$Drug.Metadata in ",
    rds_path,
    call. = FALSE
  )
}

drug_metadata <- as_plain_data_frame(mae_metadata[["Drug.Metadata"]])
drug_metadata <- collapse_list_columns(drug_metadata)
drug_metadata <- drug_metadata |>
  rownames_to_column("RDS.RowName")

dir.create(dirname(metadata_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(drug_metadata, metadata_path, na = "NA")

summary <- tibble(
  dataset = dataset,
  rows = nrow(drug_metadata),
  columns = ncol(drug_metadata),
  source_rds = rds_path,
  output = metadata_path
)
write_tsv(summary, summary_path, na = "NA")

cat(
  "[extract_sub_dataset_drug_metadata]",
  "dataset=",
  dataset,
  " rows=",
  nrow(drug_metadata),
  " columns=",
  ncol(drug_metadata),
  " output=",
  metadata_path,
  "\n",
  sep = ""
)
