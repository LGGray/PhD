library(slingshot)
library(Seurat)
library(SingleCellExperiment)

pbmc <- readRDS('pbmc.female.control-managed.RDS')

# Convert Seurat to SingleCellExperiment
sce <- as.SingleCellExperiment(pbmc)

# Slingshot
sce <- slingshot(sce, clusterLabels = 'cellTypist', reducedDim = 'UMAP')

saveRDS(sce, 'slingshot.RDS')

# Plot the trajectory
library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

pdf('slingshot.pdf')
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')
dev.off()
