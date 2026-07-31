suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(Matrix)
  library(SummarizedExperiment)
})

if (exists("snakemake")) {
  rds_path <- snakemake@input[["mae"]]
  out_dir <- snakemake@output[["outdir"]]
} else {
  args <- commandArgs(trailingOnly = TRUE)
  rds_path <- if (length(args) >= 1) args[1] else "data/results/HDD_v2.3.RDS"
  out_dir <- if (length(args) >= 2) args[2] else "data/results/HDD_v2.3_csv"
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
assay_dir <- file.path(out_dir, "assays")
dir.create(assay_dir, recursive = TRUE, showWarnings = FALSE)

mae <- readRDS(rds_path)

coldata_df <- as.data.frame(colData(mae))
write.csv(coldata_df, file.path(out_dir, "colData.csv"), row.names = TRUE)

format_sparse_values <- function(values) {
  if (is.integer(values)) {
    return(as.character(values))
  }
  format(values, scientific = FALSE, trim = TRUE)
}

write_sparse_matrix_csv <- function(mat, outfile) {
  mat <- as(mat, "RsparseMatrix")
  row_names <- rownames(mat)
  col_names <- colnames(mat)
  if (is.null(row_names)) {
    row_names <- as.character(seq_len(nrow(mat)))
  }
  if (is.null(col_names)) {
    col_names <- as.character(seq_len(ncol(mat)))
  }

  con <- file(outfile, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste(c("", col_names), collapse = ","), con = con)

  for (row_idx in seq_len(nrow(mat))) {
    row_values <- rep.int("0", ncol(mat))
    start <- mat@p[[row_idx]] + 1L
    end <- mat@p[[row_idx + 1L]]
    if (start <= end) {
      row_values[mat@j[start:end] + 1L] <- format_sparse_values(mat@x[
        start:end
      ])
    }
    writeLines(
      paste(c(row_names[[row_idx]], row_values), collapse = ","),
      con = con
    )
  }
}

exp_names <- names(experiments(mae))
for (exp_name in exp_names) {
  se <- experiments(mae)[[exp_name]]
  assay_obj <- assays(se)[[1]]
  outfile <- file.path(assay_dir, paste0(exp_name, ".csv"))
  if (inherits(assay_obj, "sparseMatrix")) {
    message("[export_mae_csvs] writing sparse assay CSV: ", outfile)
    write_sparse_matrix_csv(assay_obj, outfile)
    next
  } else if (is.data.frame(assay_obj)) {
    assay_mat <- assay_obj
  } else {
    assay_mat <- as.matrix(assay_obj)
  }

  if (is.null(rownames(assay_mat))) {
    rownames(assay_mat) <- seq_len(nrow(assay_mat))
  }
  if (is.null(colnames(assay_mat))) {
    colnames(assay_mat) <- seq_len(ncol(assay_mat))
  }

  message("[export_mae_csvs] writing assay CSV: ", outfile)
  write.csv(assay_mat, outfile, row.names = TRUE)
}
