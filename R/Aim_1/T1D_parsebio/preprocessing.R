library(Seurat)
library(SeuratDisk)
library(magrittr)
library(ddqcR)
library(DoubletFinder)

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrY.Rdata')
source('/directflow/SCCGGroupShare/projects/lacgra/R_code/functions/calc.min.pc.R')

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/T1D_parsebio')

# Convert("T1D.h5ad", dest="h5seurat", overwrite=F)
pbmc <- LoadH5Seurat("T1D.h5seurat")
metadata <- read.csv('T1D.metadata.csv', row.names=1)
colnames(metadata)[1] <- 'individual'
metadata$condition <- gsub('_[0-9]+', '', metadata$individual)
expr <- GetAssayData(pbmc, assay='RNA')
pbmc <- CreateSeuratObject(expr)
pbmc@meta.data <- cbind(pbmc@meta.data, metadata)

# Perform the filtering * Based on parsebio website *
pbmc$percent.mt <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA < 5000 & nCount_RNA < 20000 & percent.mt < 15)

# SCTransform
pbmc <- SCTransform(pbmc, verbose = FALSE)

# Export .h5ad file for cellTypist
SaveH5Seurat(pbmc, filename = "pbmc.h5Seurat", overwrite = TRUE)
Convert("pbmc.h5Seurat", dest = "h5ad", overwrite = TRUE)

mtx <- as.matrix(GetAssayData(pbmc))
write.csv(mtx, 'raw.counts.csv')

# Save RDS file for downstream cellTypist analysis
saveRDS(pbmc, 'pbmc.unlabelled.RDS')