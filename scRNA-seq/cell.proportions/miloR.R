library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(Seurat)

source('~/R_code/functions/calc.min.pc.R')

args <- commandArgs()

setwd(paste0('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/', args[1]))

pbmc <- readRDS('pbmc.female.RDS')
pbmc <- RunPCA(pbmc)
min.pc <- calc.min.pc(pbmc)
pbmc <- RunUMAP(pbmc, dims=1:min.pc)
pbmc_sce <- as.SingleCellExperiment(pbmc)
pbmc_milo <- Milo(pbmc_sce)

# Change working directory
if(dir.exists('miloR')!=T){dir.create('miloR')}
setwd('miloR')

# Construct KNN graph
traj_milo <- buildGraph(pbmc_milo, k = min.pc, d = min.pc)

# Define representative neighbourhoods
traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = min.pc, d=min.pc, refined = TRUE)

pdf('NhoodSizeHist.pdf')
plotNhoodSizeHist(traj_milo)
dev.off()

# Counting cells in neighbourhoods
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), 
                        sample="individual")

# Differential abundance testing
traj_design <- data.frame(colData(traj_milo))[,c("individual", "condition")]
traj_design$condition <- factor(traj_design$condition, levels = c('SLE', 'Control'))
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$individual

# Spatial FDR correction. Store distances between nearest neighbors
traj_milo <- calcNhoodDistance(traj_milo, d=min.pc)

# Differential abundance testing
da_results <- testNhoods(traj_milo, design = ~ condition, design.df = traj_design)

da_results %>%
  arrange(SpatialFDR) %>%
  head()

# DA testing
da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "predicted.celltype.l2")
# Mark mixed celltypes
da_results$predicted.celltype.l2 <- ifelse(da_results$predicted.celltype.l2_fraction < 0.7, "Mixed", da_results$predicted.celltype.l2)

da_results <- groupNhoods(traj_milo, da_results, max.lfc.delta = 2)

save(da_results, file='da_results.Rdata')

# Nhood graph
traj_milo <- buildNhoodGraph(traj_milo)
pdf('NhoodGroups.pdf')
plotNhoodGroups(traj_milo, da_results, layout="UMAP")
dev.off()

pdf('pvalue.hist.pdf')
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
dev.off()

pdf('volcano.plot.pdf')
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1)
dev.off()

# Distribution of DA fold changes
pdf('DA.beeswarm.pdf')
plotDAbeeswarm(da_results, group.by = "predicted.celltype.l2")
dev.off()

# Visualise neighbourhoods displaying DA
pdf('NhoodGraph.pdf')
traj_milo <- buildNhoodGraph(traj_milo)
dev.off()

# cellcount <- pbmc@meta.data %>% 
#   group_by(condition) %>%
#   count(predicted.celltype.l2)

saveRDS(traj_milo, 'traj_milo.RDS')