library(MultiAssayExperiment)
library(Matrix)
library(SummarizedExperiment)


colData <- read.csv(
  "data/procdata/colData.csv",
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
rownames(colData) <- colData$Pubchem.CID
colData <- DataFrame(colData, row.names = rownames(colData))

bioassays <- read.csv(
  "data/procdata/experiments/bioassays.csv",
  row.names = 1,
  check.names = FALSE,
  na.strings = c("Not Measured")
)
toxcast <- read.csv(
  "data/procdata/experiments/toxcast.csv",
  row.names = 1,
  check.names = FALSE
)
colnames(toxcast) <- sub("\\.0$", "", as.character(colnames(toxcast)))

tox21 <- read.csv(
  "data/procdata/experiments/tox21.csv",
  row.names = 1,
  check.names = FALSE
)
clintox <- read.csv(
  "data/procdata/experiments/clintox.csv",
  row.names = 1,
  check.names = FALSE
)
sider <- read.csv(
  "data/procdata/experiments/sider.csv",
  row.names = 1,
  check.names = FALSE
)


fingerprint.files <- list.files("data/procdata/experiments/fingerprints/")
fp_assays <- list()

for (fp.file in fingerprint.files) {
  # print(paste0("data/procdata/experiments/fingerprints/",fp.file))
  fp.data <- read.csv(
    paste0("data/procdata/experiments/fingerprints/", fp.file),
    check.names = FALSE
  )
  fp.data <- as(fp.data, "sparseMatrix")

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
saveRDS(mae, "data/results/HDD_v1.RDS")
print(mae)
