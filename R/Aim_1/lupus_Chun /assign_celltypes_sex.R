# Script for building Seurat object from UC_GSE125527 data

library(Seurat)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
args <- as.numeric(args)

ancestry <- c('european', 'asian')

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrY.Rdata')

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/lupus_Chun')

pbmc <- readRDS(paste0('pbmc.', ancestry[args], '.unlabelled.RDS'))

labels <- read.csv(paste0('cellTypist/predicted_labels.', ancestry[args], '.csv'))
pbmc@meta.data <- cbind(pbmc@meta.data, cellTypist=labels$majority_voting)

Idents(pbmc) <- 'cellTypist'

pdf(paste('DimPlot.cellTypist.female.', ancestry[args], '.pdf'))
DimPlot(pbmc, label = TRUE, reduction='umap', repel=T)
dev.off()

saveRDS(pbmc, paste0('pbmc.female', ancestry[args], '.RDS'))