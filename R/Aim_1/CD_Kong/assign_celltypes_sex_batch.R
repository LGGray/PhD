# Script for building Seurat object from UC_GSE182270 data

library(Seurat)
library(magrittr)

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrY.Rdata')

pbmc <- readRDS('pbmc.unlabelled.RDS')

labels <- read.csv('cellTypist/predicted_labels.csv')
pbmc@meta.data <- cbind(pbmc@meta.data, cellTypist=labels$majority_voting)

Idents(pbmc) <- 'cellTypist'

# Load in IA-SVA output
iasva <- readRDS("iasva.res.RDS")
pbmc$SV1 <- iasva$sv[,1]
pbmc$SV2 <- iasva$sv[,2]

# Add sex to metadata
pbmc$sex <- ifelse(pbmc$sex == 'male', 'M', 'F')

# Output all cells
saveRDS(pbmc, 'pbmc.RDS')

# Subset for females and output 
pbmc.female <- subset(pbmc, sex == 'F')

pdf('DimPlot.female.pdf')
DimPlot(pbmc.female, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()
saveRDS(pbmc.female, 'pbmc.female.RDS')

# Subset for males and output 
pbmc.male <- subset(pbmc, sex == 'M')

pdf('DimPlot.male.pdf')
DimPlot(pbmc.male, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()
saveRDS(pbmc.male, 'pbmc.male.RDS')

