# Script for building Seurat object from lupus_Chun data
library(Seurat)
library(magrittr)

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/lupus_Chun')

pbmc <- readRDS('pbmc.unlabelled.RDS')

labels <- read.csv('cellTypist/predicted_labels.csv')
pbmc@meta.data <- cbind(pbmc@meta.data, cellTypist=labels$majority_voting)

Idents(pbmc) <- 'cellTypist'

# Load in IA-SVA output
iasva <- readRDS("iasva.res.RDS")
pbmc$SV1 <- iasva$sv[,1]

# Output all cells
pdf('DimPlot.cellTypist.all.pdf')
DimPlot(pbmc, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()
saveRDS(pbmc, 'pbmc.RDS')

# Split by sex and save plots and objects
pbmc.male <- subset(pbmc, sex == 'M')
pdf('DimPlot.cellTypist.male.pdf')
DimPlot(pbmc.male, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()
saveRDS(pbmc.male, 'pbmc.male.RDS')

pbmc.female <- subset(pbmc, sex == 'F')
pdf('DimPlot.cellTypist.female.pdf')
DimPlot(pbmc.female, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()
saveRDS(pbmc, 'pbmc.female.RDS')