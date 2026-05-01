suppressPackageStartupMessages({
  library(data.table)
  library(MultiAssayExperiment)
  library(Matrix)
  library(SummarizedExperiment)
})

read_table <- function(path, na.strings = c("NA", "", "None", "Unknown", "-")) {
  data.table::fread(
    path,
    check.names = FALSE,
    data.table = FALSE,
    na.strings = na.strings
  )
}

read_assay_table <- function(path, na.strings = c("NA", "")) {
  assay <- read_table(path, na.strings = na.strings)
  row_ids <- assay[[1]]
  assay[[1]] <- NULL
  rownames(assay) <- row_ids
  assay
}

normalize_logical_column <- function(df, col_name) {
  if (!(col_name %in% colnames(df))) {
    return(df)
  }

  values <- df[[col_name]]
  if (is.logical(values)) {
    return(df)
  }

  normalized <- toupper(trimws(as.character(values)))
  if (!all(is.na(values) | normalized %in% c("TRUE", "FALSE"))) {
    return(df)
  }

  df[[col_name]] <- normalized == "TRUE"
  df[[col_name]][is.na(values)] <- NA
  df
}

if (exists("snakemake")) {
  coldata_path <- snakemake@input[["colData"]]
  bioassays_path <- snakemake@input[["bioassays"]]
  toxcast_path <- snakemake@input[["toxcast"]]
  tox21_path <- snakemake@input[["tox21"]]
  clintox_path <- snakemake@input[["clintox"]]
  sider_path <- snakemake@input[["sider"]]
  fingerprint_files <- as.character(snakemake@input[["fingerprints"]])
  output_path <- snakemake@output[["mae"]]
} else {
  coldata_path <- "data/procdata/colData.csv"
  bioassays_path <- "data/procdata/experiments/bioassays.csv"
  toxcast_path <- "data/procdata/experiments/toxcast.csv"
  tox21_path <- "data/procdata/experiments/tox21.csv"
  clintox_path <- "data/procdata/experiments/clintox.csv"
  sider_path <- "data/procdata/experiments/sider.csv"
  fingerprint_files <- list.files(
    "data/procdata/experiments/fingerprints/",
    full.names = TRUE
  )
  output_path <- "data/results/HDD_v2.RDS"
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

colData <- read_table(coldata_path)
for (logical_col in c(
  "In.AnnotationDB",
  "FDA.Approved",
  "In.LINCS",
  "In.JUMP.CP",
  "In.OASIS",
  "In.GEOM"
)) {
  colData <- normalize_logical_column(colData, logical_col)
}
if (!("HDD.Compound.ID" %in% colnames(colData))) {
  stop("colData is missing HDD.Compound.ID", call. = FALSE)
}
if (
  any(is.na(colData$HDD.Compound.ID)) || anyDuplicated(colData$HDD.Compound.ID)
) {
  stop("HDD.Compound.ID must be unique and non-missing", call. = FALSE)
}
rownames(colData) <- colData$HDD.Compound.ID
colData <- DataFrame(colData, row.names = rownames(colData))

bioassays <- read_assay_table(
  bioassays_path,
  na.strings = c("Not Measured")
)
toxcast <- read_assay_table(toxcast_path)

tox21 <- read_assay_table(tox21_path)
clintox <- read_assay_table(clintox_path)
sider <- read_assay_table(sider_path)


fp_assays <- list()

for (fingerprint_file in fingerprint_files) {
  fp.file <- basename(fingerprint_file)
  if (grepl("\\.mtx$", fp.file)) {
    fp.data <- readMM(fingerprint_file)
    fp.data <- as(fp.data, "CsparseMatrix")
    rownames(fp.data) <- paste0("V", seq_len(nrow(fp.data)))
    colnames(fp.data) <- rownames(colData)
  } else {
    fp.data <- read.csv(
      fingerprint_file,
      check.names = FALSE
    )
    fp.data <- as(fp.data, "sparseMatrix")
  }

  fp.stem <- unlist(strsplit(fp.file, "\\."))
  fp.name <- paste(
    c("fingerprint", as.character(fp.stem[-length(fp.stem)])),
    sep = "",
    collapse = "."
  )
  fp_assays[[fp.name]] <- SummarizedExperiment(assays = list(fp.name = fp.data))
}


experiments <- c(
  list(
    SIDER = SummarizedExperiment(assays = list(SIDER = as.matrix(sider))),
    Bioassays = SummarizedExperiment(assays = list(Bioassays = bioassays)),
    Tox21 = SummarizedExperiment(assays = list(Tox21 = tox21)),
    ToxCast = SummarizedExperiment(assays = list(ToxCast = toxcast)),
    ClinTox = SummarizedExperiment(assays = list(ClinTox = clintox))
  ),
  fp_assays
)


experimentList <- as(experiments, "ExperimentList")

sampleMapList <- lapply(experiments, function(se) {
  data.frame(
    primary = colnames(se),
    colname = colnames(se),
    stringsAsFactors = FALSE
  )
})

mae <- MultiAssayExperiment(
  experiments = experimentList,
  colData = colData,
  sampleMap = listToMap(sampleMapList)
)
saveRDS(mae, output_path)
print(mae)
