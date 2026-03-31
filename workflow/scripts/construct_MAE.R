suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(Matrix)
  library(SummarizedExperiment)
})

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
  fingerprint_files <- list.files("data/procdata/experiments/fingerprints/", full.names = TRUE)
  output_path <- "data/results/HDD_v1.1.RDS"
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

colData <- read.csv(
  coldata_path,
  na.strings = c("NA", "", "None", "Unknown", "-")
)
colnames(colData) <- sub(
  "^Hepatotoxicity\\.Likelihood\\.\\.Detailed\\.$",
  "Hepatotoxicity.Likelihood.Detailed",
  colnames(colData)
)
colnames(colData) <- sub(
  "^Hepatotoxiciy\\.Likelihood\\.\\.Score\\.$",
  "Hepatotoxicity.Likelihood.Score",
  colnames(colData)
)
for (logical_col in c("FDA.Approved", "In.L1000", "In.JUMP.CP")) {
  colData <- normalize_logical_column(colData, logical_col)
}
rownames(colData) <- colData$Pubchem.CID
colData <- DataFrame(colData, row.names = rownames(colData))

bioassays <- read.csv(
  bioassays_path,
  row.names = 1,
  check.names = FALSE,
  na.strings = c("Not Measured")
)
toxcast <- read.csv(
  toxcast_path,
  row.names = 1,
  check.names = FALSE
)
colnames(toxcast) <- sub("\\.0$", "", as.character(colnames(toxcast)))

tox21 <- read.csv(
  tox21_path,
  row.names = 1,
  check.names = FALSE
)
clintox <- read.csv(
  clintox_path,
  row.names = 1,
  check.names = FALSE
)
sider <- read.csv(
  sider_path,
  row.names = 1,
  check.names = FALSE
)


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

mae = MultiAssayExperiment(
  experiments = experimentList,
  colData = colData,
  sampleMap = listToMap(sampleMapList)
)
saveRDS(mae, output_path)
print(mae)
