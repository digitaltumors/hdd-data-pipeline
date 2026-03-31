suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(Matrix)
})

if (exists("snakemake")) {
  rds_path <- snakemake@input[["mae"]]
  out_dir <- snakemake@output[["outdir"]]
} else {
  args <- commandArgs(trailingOnly = TRUE)
  rds_path <- if (length(args) >= 1) args[1] else "data/results/HDD_v1.1.RDS"
  out_dir <- if (length(args) >= 2) args[2] else "data/results/HDD_v1.1_csv"
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
assay_dir <- file.path(out_dir, "assays")
dir.create(assay_dir, recursive = TRUE, showWarnings = FALSE)

mae <- readRDS(rds_path)

coldata_df <- as.data.frame(colData(mae))
write.csv(coldata_df, file.path(out_dir, "colData.csv"), row.names = TRUE)

exp_names <- names(experiments(mae))
for (exp_name in exp_names) {
  se <- experiments(mae)[[exp_name]]
  assay_obj <- assays(se)[[1]]
  if (inherits(assay_obj, "sparseMatrix")) {
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

  outfile <- file.path(assay_dir, paste0(exp_name, ".csv"))
  write.csv(assay_mat, outfile, row.names = TRUE)
}
