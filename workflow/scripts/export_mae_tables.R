suppressPackageStartupMessages({
  library(Matrix)
  library(MultiAssayExperiment)
  library(S4Vectors)
  library(SummarizedExperiment)
})

write_tsv <- function(frame, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    frame,
    path,
    sep = "\t",
    quote = TRUE,
    qmethod = "double",
    row.names = FALSE,
    na = "NA"
  )
}

flatten_config <- function(value, path = "") {
  if (is.list(value)) {
    item_names <- names(value)
    if (is.null(item_names)) {
      item_names <- sprintf("[%d]", seq_along(value))
    }
    rows <- lapply(seq_along(value), function(index) {
      separator <- if (startsWith(item_names[[index]], "[")) "" else "."
      item_path <- if (nzchar(path)) {
        paste0(path, separator, item_names[[index]])
      } else {
        item_names[[index]]
      }
      flatten_config(value[[index]], item_path)
    })
    return(do.call(rbind, rows))
  }

  value_length <- length(value)
  data.frame(
    Key = if (value_length > 1L) {
      paste0(path, "[", seq_len(value_length), "]")
    } else {
      path
    },
    Value = if (value_length == 0L) NA_character_ else as.character(value),
    stringsAsFactors = FALSE
  )
}

if (exists("snakemake")) {
  rds_path <- snakemake@input[["mae"]]
  out_dir <- snakemake@output[["outdir"]]
} else {
  rds_path <- "data/results/HDD_v3.RDS"
  out_dir <- "data/results/HDD_v3_tables"
}

if (dir.exists(out_dir)) {
  unlink(out_dir, recursive = TRUE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mae <- readRDS(rds_path)
inventory <- list()

record_file <- function(path, category, rows, columns) {
  inventory[[length(inventory) + 1L]] <<- data.frame(
    File = sub(paste0("^", out_dir, "/?"), "", path),
    Category = category,
    Rows = rows,
    Columns = columns,
    Bytes = file.info(path)$size,
    SHA256 = unname(tools::sha256sum(path)),
    stringsAsFactors = FALSE
  )
}

object_metadata <- S4Vectors::metadata(mae)
pipeline_config <- flatten_config(object_metadata$Pipeline$Config)
pipeline_config_path <- file.path(
  out_dir,
  "metadata",
  "pipeline_config.tsv"
)
write_tsv(pipeline_config, pipeline_config_path)
record_file(
  pipeline_config_path,
  "config",
  nrow(pipeline_config),
  ncol(pipeline_config)
)

coldata <- as.data.frame(MultiAssayExperiment::colData(mae))
coldata_path <- file.path(out_dir, "colData.tsv")
write_tsv(coldata, coldata_path)
record_file(coldata_path, "colData", nrow(coldata), ncol(coldata))

indications <- as.data.frame(object_metadata$Drug.Indications)
indications_path <- file.path(out_dir, "metadata", "drug_indications.tsv")
write_tsv(indications, indications_path)
record_file(
  indications_path,
  "metadata",
  nrow(indications),
  ncol(indications)
)

sample_map <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
sample_map_path <- file.path(out_dir, "sampleMap.tsv")
write_tsv(sample_map, sample_map_path)
record_file(sample_map_path, "sampleMap", nrow(sample_map), ncol(sample_map))

for (experiment_name in names(MultiAssayExperiment::experiments(mae))) {
  experiment <- MultiAssayExperiment::experiments(mae)[[experiment_name]]
  assay_value <- SummarizedExperiment::assay(experiment, 1L)
  assay_dir <- file.path(out_dir, "assays", experiment_name)
  dir.create(assay_dir, recursive = TRUE, showWarnings = FALSE)

  if (inherits(assay_value, "sparseMatrix")) {
    matrix_path <- file.path(assay_dir, "matrix.mtx")
    Matrix::writeMM(assay_value, matrix_path)
    record_file(
      matrix_path,
      "sparse_assay",
      nrow(assay_value),
      ncol(assay_value)
    )

    rows_path <- file.path(assay_dir, "rows.tsv")
    columns_path <- file.path(assay_dir, "columns.tsv")
    write_tsv(
      data.frame(
        Row = seq_len(nrow(assay_value)),
        Feature = rownames(assay_value)
      ),
      rows_path
    )
    write_tsv(
      data.frame(
        Column = seq_len(ncol(assay_value)),
        HDD.Compound.ID = colnames(assay_value)
      ),
      columns_path
    )
    record_file(rows_path, "sparse_row_map", nrow(assay_value), 2L)
    record_file(columns_path, "sparse_column_map", ncol(assay_value), 2L)
    next
  }

  dense_path <- file.path(assay_dir, "matrix.tsv")
  dense_frame <- data.frame(
    Feature = rownames(assay_value),
    as.data.frame(assay_value, check.names = FALSE),
    check.names = FALSE
  )
  write_tsv(dense_frame, dense_path)
  record_file(
    dense_path,
    "dense_assay",
    nrow(assay_value),
    ncol(assay_value)
  )
}

manifest <- do.call(rbind, inventory)
write_tsv(manifest, file.path(out_dir, "file_manifest.tsv"))
