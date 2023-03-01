# Script for building Seurat object from lupus_James data
library(Seurat)
library(magrittr)

ancestry = commandArgs(trailingOnly=TRUE)

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/lupus_Chun')

pbmc <- readRDS(paste0('pbmc.', ancestry, '.unlabelled.RDS'))

labels <- read.csv(paste0('cellTypist.', ancestry, '/predicted_labels.csv'))
pbmc@meta.data <- cbind(pbmc@meta.data, cellTypist=labels$majority_voting)

Idents(pbmc) <- 'cellTypist'

# split data into male and female
pbmc.female <- subset(pbmc, sex =='female')

pdf(paste('DimPlot.female.', ancestry, '.pdf'))
DimPlot(pbmc.female, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()
saveRDS(pbmc.female, paste0('pbmc.female.', ancestry, '.RDS'))

pbmc.male <- subset(pbmc, sex == 'male')

pdf(paste('DimPlot.male.', ancestry, '.pdf'))
DimPlot(pbmc.male, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()
saveRDS(pbmc.male, paste0('pbmc.male.', ancestry, '.RDS'))