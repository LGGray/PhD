# Script for building Seurat object from UC_GSE182270 data

library(Seurat)
library(magrittr)

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/UC_JamesLab')

pbmc <- readRDS('pbmc.unlabelled.RDS')

labels <- read.csv('cellTypist/predicted_labels.csv')
pbmc@meta.data <- cbind(pbmc@meta.data, cellTypist=labels$majority_voting)

Idents(pbmc) <- 'cellTypist'

pdf('DimPlot.cellTypist.all.pdf')
DimPlot(pbmc, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()

# Output all cells
saveRDS(pbmc, 'pbmc.RDS')

# Subset for females and output 
pbmc.female <- subset(pbmc, sex == 'female')

pdf('DimPlot.female.pdf')
DimPlot(pbmc.female, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()

saveRDS(pbmc.female, 'pbmc.female.RDS')

# Subset for males and output 
pbmc.male <- subset(pbmc, sex == 'male')

pdf('DimPlot.male.pdf')
DimPlot(pbmc.male, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()

saveRDS(pbmc.male, 'pbmc.male.RDS')