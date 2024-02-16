library(slingshot)
library(Seurat)
library(SingleCellExperiment)
library(tradeSeq)

pbmc <- readRDS('pbmc.female.control-managed.RDS')

T_cells <- subset(pbmc, idents = c("Tcm/Naive helper T cells", "Tem/Effector helper T cells",
"Tem/Temra cytotoxic T cells", "Tem/Trm cytotoxic T cells", "Tcm/Naive cytotoxic T cells",
"MAIT cells", "Cycling T cells", "Regulatory T cells"))

# Convert Seurat to SingleCellExperiment
sce <- as.SingleCellExperiment(T_cells)

# Slingshot
sce <- slingshot(sce, clusterLabels = 'cellTypist', reducedDim = 'UMAP', approx_points = 150, omega = TRUE)

saveRDS(sce, 'slingshot_T.cells.RDS')

# Plot the trajectory
library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

pdf('slingshot.pdf')
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
dev.off()
