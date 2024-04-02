library(celda)
library(SingleCellExperiment)
library(Seurat)

pbmc <- readRDS('pbmc.female.RDS')
pbmc@assays$RNA@scale.data <- as.matrix(0)

# convert to SingleCellExperiment
sce <- as.SingleCellExperiment(pbmc)

sce <- selectFeatures(sce)

moduleSplit <- recursiveSplitModule(sce, initialL = 2, maxL = 15)

pdf('perplexity.elbow.pdf')
plotGridSearchPerplexity(moduleSplit)
dev.off()

pdf('PRC.pdf')
plotRPC(moduleSplit)
dev.off()

sce <- celda_CG(x = sce, K = length(levels(pbmc)), L = 10, verbose = FALSE, nchains = 3)

pdf('celdaHeatmap.pdf')
plot(celdaHeatmap(sce = sce, nfeatures = 10))
dev.off()

geneSetEnrich(sce,
  databases = c("GO_Biological_Process_2018"))


load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
modules <- celdaModules(sce)
names(modules) <- rownames(altExp(sce))

modules[names(modules) %in% rownames(chrX)]

names(modules[modules == 2])