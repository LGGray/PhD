library(Seurat)
library(magrittr)
library(ddqcR)
library(reticulate)
library(celda)
library(BiocParallel)
library(scDblFinder)
library(transformGamPoi)
library(iasva)
library(sva)
library(irlba)
library(SummarizedExperiment)

options(future.globals.maxSize = 891289600)

# Read in Seurat object
pbmc <- readRDS('pbmc.female.RDS')
DefaultAssay(pbmc) <- 'RNA'

# Subset for flare and managed samples
pbmc.disease_state <- subset(pbmc, disease_state %in% c('flare', 'managed'))

# Normalise data with Delta method-based variance stabilizing
exp.matrix <- GetAssayData(pbmc.disease_state, slot = 'counts')
exp.matrix.transformed <- acosh_transform(exp.matrix)

# Add transformed data to Seurat object
pbmc.disease_state <- SetAssayData(object=pbmc.disease_state, slot = 'counts', new.data=exp.matrix.transformed)

# Cell type clustering and inspection of markers
pbmc.disease_state <- FindVariableFeatures(pbmc.disease_state)
pbmc.disease_state <- ScaleData(pbmc.disease_state)
pbmc.disease_state <- RunPCA(pbmc.disease_state)

# Determine percent of variation associated with each PC
pct <- pbmc[["pca"]]@stdev / sum(pbmc[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# Minimum of the two calculation
pcs <- min(co1, co2)
print(paste('Selected # PCs', pcs))

pbmc <- FindNeighbors(pbmc, dims=1:pcs)
pbmc <- FindClusters(pbmc, resolution=1)
pbmc <- RunUMAP(pbmc, dims = 1:pcs)

Idents(pbmc.disease_state) <- 'cellTypist'

pdf('disease_state/DimPlot.cellTypist.disease_state.pdf')
DimPlot(pbmc.disease_state, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()

saveRDS(pbmc.disease_state, 'disease_state/pbmc.female.RDS')