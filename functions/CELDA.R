library(celda)
library(SingleCellExperiment)
library(Seurat)

pbmc <- readRDS('pbmc.female.control-managed.RDS')

ncM <- subset(pbmc, idents = c('Non-classical monocytes'))

# convert to SingleCellExperiment
sce <- as.SingleCellExperiment(ncM)

sce <- selectFeatures(sce)

# identify best K and L
moduleSplit <- recursiveSplitModule(sce, initialL = 2, maxL = 15)

pdf('moduleSplit.pdf')
plotGridSearchPerplexity(moduleSplit)
dev.off()

# sce <- celda_CG(x = simsce, K = , L = 10, verbose = FALSE, nchains = 1)

# # Heatmap
# plot(celdaHeatmap(sce = sce, nfeatures = 10))

# # relationship between modules and cell types
# celdaProbabilityMap(sce)
# dev.off()

# # Coexpression heatmap
# moduleHeatmap(sce, featureModule = c(1,2), topCells = 100)
# dev.off()