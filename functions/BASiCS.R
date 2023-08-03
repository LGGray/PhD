library(BASiCS)?
# library(Seurat)

pbmc <- readRDS('pbmc.female.RDS')
control <- subset(pbmc, cellTypist == 'Regulatory T cells' & condition == 'control')
control <- subset(control, features=names(rowSums(control) > 0))
control_sce <- SingleCellExperiment(assays = list(counts = control@assays$RNA@counts))
colData(control_sce)$BatchInfo <- control$SV1

disease <- subset(pbmc, cellTypist == 'Regulatory T cells' & condition == 'disease')
disease <- subset(disease, features=names(rowSums(disease) > 0))
disease_sce <- SingleCellExperiment(assays = list(counts = disease@assays$RNA@counts))
colData(disease_sce)$BatchInfo <- disease$SV1

set.seed(42)
control_Chain <- BASiCS_MCMC(
  Data = control_sce,
  N = 20000, Thin = 20, Burn = 10000,
  PrintProgress = FALSE, Regression = TRUE,
  WithSpikes = FALSE,
  StoreChains = TRUE, StoreDir = 'BASiCS', RunName = "Regulatory.T.cells.control"
)
disease_Chain <- BASiCS_MCMC(
  Data = disease_sce,
  N = 20000, Thin = 20, Burn = 1000,
  PrintProgress = FALSE, Regression = TRUE,
  WithSpikes = FALSE,
  StoreChains = TRUE, StoreDir = 'BASiCS', RunName = "Regulatory.T.cells.disease"
)

# If Acceptance rates are above 0.44 increase N and Burn

# Convergence assessment

# Differential expression
Test <- BASiCS_TestDE(
  Chain1 = control_Chain, Chain2 = disease_Chain,
  GroupLabel1 = "control", GroupLabel2 = "disease",
  EpsilonM = log2(1.5), EpsilonD = log2(1.5),
  EFDR_M = 0.10, EFDR_D = 0.10,
  Offset = TRUE, PlotOffset = TRUE, Plot = TRUE
)
saveRDS(Test,"BASiCS/Regulatory.T.cells_TestDE.RDS")

BASiCS_PlotDE(Test)