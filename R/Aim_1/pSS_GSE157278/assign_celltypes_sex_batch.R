# Script for building Seurat object from UC_GSE182270 data

library(Seurat)
library(magrittr)

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrY.Rdata')

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/pSS_GSE157278')

pbmc <- readRDS('pbmc.unlabelled.RDS')

labels <- read.csv('cellTypist/predicted_labels.csv')
pbmc@meta.data <- cbind(pbmc@meta.data, cellTypist=labels$majority_voting)

Idents(pbmc) <- 'cellTypist'

# Load in SVA output
load("svaseq.RData")
pbmc$batch_1 <- svseq$sv[,1]
pbmc$batch_2 <- svseq$sv[,2]

# Add sex to metadata
pbmc$sex <- 'F'

pdf('DimPlot.cellTypist.female.pdf')
DimPlot(pbmc, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()

# Output all cells
saveRDS(pbmc, 'pbmc.female.RDS')
