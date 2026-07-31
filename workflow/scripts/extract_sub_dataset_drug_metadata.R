suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(readr)
  library(S4Vectors)
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

sanitize_control_characters <- function(df) {
  for (column in names(df)) {
    if (!is.character(df[[column]]) && !is.factor(df[[column]])) {
      next
    }

    values <- as.character(df[[column]])
    values <- gsub("[\t\r\n]+", " ", values)
    df[[column]] <- trimws(values)
  }
  df
}

first_non_missing_column <- function(df, candidates) {
  values <- rep(NA_character_, nrow(df))
  for (column in candidates) {
    if (!(column %in% colnames(df))) {
      next
    }

    candidate <- trimws(as.character(df[[column]]))
    candidate[
      is.na(candidate) |
        candidate == "" |
        toupper(candidate) %in% c("NA", "NAN", "NONE", "NULL")
    ] <- NA_character_
    use <- is.na(values) & !is.na(candidate)
    values[use] <- candidate[use]
  }
  values
}

normalize_identity_columns <- function(
  df,
  source_key_column,
  source_key_columns
) {
  if (length(source_key_columns) > 0) {
    df[[source_key_column]] <- first_non_missing_column(
      df,
      c(source_key_column, source_key_columns)
    )
  }

  df[["Pubchem.CID"]] <- first_non_missing_column(
    df,
    c(
      "Pubchem.CID",
      "PubChem.CID",
      "AnnotationDB.PubChem.CID",
      "AnnotationGx.PubChem.CID",
      "cid"
    )
  )
  df[["InChIKey"]] <- first_non_missing_column(
    df,
    c(
      "InChIKey",
      "AnnotationDB.InChIKey",
      "AnnotationGx.InChIKey",
      "inchikey"
    )
  )
  df[["SMILES"]] <- first_non_missing_column(
    df,
    c(
      "SMILES",
      "AnnotationDB.Canonical.SMILES",
      "AnnotationDB.Connectivity.SMILES",
      "AnnotationGx.Connectivity.SMILES",
      "AnnotationGx.Canonical.SMILES",
      "cpd_smiles",
      "smiles"
    )
  )
  df
}

if (exists("snakemake")) {
  dataset <- snakemake@params[["dataset"]]
  object_type <- snakemake@params[["object_type"]]
  source_key_column <- snakemake@params[["source_key_column"]]
  source_key_columns <- as.character(
    unlist(snakemake@params[["source_key_source_columns"]])
  )
  rds_path <- snakemake@input[["rds"]]
  metadata_path <- snakemake@output[["metadata"]]
  summary_path <- snakemake@output[["summary"]]
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 7) {
    stop(
      "Usage: extract_sub_dataset_drug_metadata.R ",
      "<dataset> <object_type> <source_key_column> ",
      "<source_key_columns_csv> <input_rds> <metadata_tsv> <summary_tsv>",
      call. = FALSE
    )
  }
  dataset <- args[[1]]
  object_type <- args[[2]]
  source_key_column <- args[[3]]
  source_key_columns <- strsplit(args[[4]], ",", fixed = TRUE)[[1]]
  rds_path <- args[[5]]
  metadata_path <- args[[6]]
  summary_path <- args[[7]]
}

dataset_object <- readRDS(rds_path)
if (identical(object_type, "mae")) {
  if (!inherits(dataset_object, "MultiAssayExperiment")) {
    stop("Configured MAE input is not a MultiAssayExperiment: ", rds_path)
  }
  object_metadata <- metadata(dataset_object)
  if (!("Drug.Metadata" %in% names(object_metadata))) {
    stop(
      "[extract_sub_dataset_drug_metadata] missing ",
      "metadata(mae)$Drug.Metadata in ",
      rds_path,
      call. = FALSE
    )
  }
  drug_metadata <- as_plain_data_frame(object_metadata[["Drug.Metadata"]])
} else if (identical(object_type, "pset")) {
  if (!identical(as.character(attr(dataset_object, "class")), "PharmacoSet")) {
    stop("Configured PSet input is not a PharmacoSet: ", rds_path)
  }
  drug_metadata <- as_plain_data_frame(attr(dataset_object, "treatment"))
} else {
  stop("Unsupported sub-dataset object_type: ", object_type, call. = FALSE)
}

drug_metadata <- collapse_list_columns(drug_metadata)
drug_metadata <- sanitize_control_characters(drug_metadata)
drug_metadata <- normalize_identity_columns(
  drug_metadata,
  source_key_column,
  source_key_columns
)
drug_metadata <- drug_metadata |>
  rownames_to_column("RDS.RowName")

dir.create(dirname(metadata_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(drug_metadata, metadata_path, na = "NA")

summary <- tibble(
  dataset = dataset,
  object_type = object_type,
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
