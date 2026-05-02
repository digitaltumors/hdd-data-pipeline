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
  bindingdb_path <- snakemake@input[["binding_db"]]
  toxcast_path <- snakemake@input[["toxcast"]]
  tox21_path <- snakemake@input[["tox21"]]
  clintox_path <- snakemake@input[["clintox"]]
  sider_path <- snakemake@input[["sider"]]
  fingerprint_files <- as.character(snakemake@input[["fingerprints"]])
  fingerprint_columns_path <- snakemake@input[["fingerprint_columns"]]
  output_path <- snakemake@output[["mae"]]
} else {
  coldata_path <- "data/procdata/colData.csv"
  bioassays_path <- "data/procdata/experiments/bioassays.csv"
  bindingdb_path <- "data/procdata/experiments/binding_db.csv"
  toxcast_path <- "data/procdata/experiments/toxcast.csv"
  tox21_path <- "data/procdata/experiments/tox21.csv"
  clintox_path <- "data/procdata/experiments/clintox.csv"
  sider_path <- "data/procdata/experiments/sider.csv"
  fingerprint_files <- list.files(
    "data/procdata/experiments/fingerprints/",
    pattern = "\\.mtx$",
    full.names = TRUE
  )
  fingerprint_columns_path <- "data/procdata/experiments/fingerprints/fingerprint_columns.tsv"
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

fingerprint_columns <- read_table(fingerprint_columns_path)
required_fingerprint_columns <- c("Fingerprint.Column", "HDD.Compound.ID")
missing_fingerprint_columns <- setdiff(
  required_fingerprint_columns,
  colnames(fingerprint_columns)
)
if (length(missing_fingerprint_columns) > 0) {
  stop(
    paste(
      "Fingerprint column map is missing columns:",
      paste(missing_fingerprint_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}
fingerprint_columns <- fingerprint_columns[
  order(fingerprint_columns$Fingerprint.Column),
  ,
  drop = FALSE
]
if (
  !identical(
    as.integer(fingerprint_columns$Fingerprint.Column),
    seq_len(nrow(fingerprint_columns))
  )
) {
  stop(
    "Fingerprint column map must use contiguous 1-based column indexes",
    call. = FALSE
  )
}
fingerprint_colnames <- as.character(fingerprint_columns$HDD.Compound.ID)
if (any(is.na(fingerprint_colnames)) || anyDuplicated(fingerprint_colnames)) {
  stop(
    "Fingerprint HDD.Compound.ID values must be unique and non-missing",
    call. = FALSE
  )
}
unknown_fingerprint_ids <- setdiff(fingerprint_colnames, rownames(colData))
if (length(unknown_fingerprint_ids) > 0) {
  stop("Fingerprint column map contains IDs absent from colData", call. = FALSE)
}

bioassays <- read_assay_table(
  bioassays_path,
  na.strings = c("Not Measured")
)
bindingdb <- read_assay_table(bindingdb_path)
toxcast <- read_assay_table(toxcast_path)

tox21 <- read_assay_table(tox21_path)
clintox <- read_assay_table(clintox_path)
sider <- read_assay_table(sider_path)


fp_assays <- list()

for (fingerprint_file in fingerprint_files) {
  fp.file <- basename(fingerprint_file)
  if (!grepl("\\.mtx$", fp.file)) {
    stop(
      paste("Unsupported fingerprint assay format:", fingerprint_file),
      call. = FALSE
    )
  }

  fp.data <- readMM(fingerprint_file)
  fp.data <- as(fp.data, "CsparseMatrix")
  if (ncol(fp.data) != length(fingerprint_colnames)) {
    stop(
      paste(
        "Fingerprint matrix column count does not match column map:",
        fingerprint_file
      ),
      call. = FALSE
    )
  }
  rownames(fp.data) <- paste0("V", seq_len(nrow(fp.data)))
  colnames(fp.data) <- fingerprint_colnames

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
    BindingDB = SummarizedExperiment(
      assays = list(BindingDB = as.matrix(bindingdb))
    ),
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

mae <- suppressMessages(MultiAssayExperiment(
  experiments = experimentList,
  colData = colData,
  sampleMap = listToMap(sampleMapList)
))

# The constructor harmonizes away colData rows that are absent from every assay.
# HDD keeps those compounds in colData, even when no parseable SMILES is available
# for fingerprint generation.
slot(mae, "colData") <- colData
stopifnot(validObject(mae))

saveRDS(mae, output_path)
print(mae)
