library(Seurat)

# read in chrX
load('../datasets/XCI/chrX.Rdata')

# Read in each seurat object
AD <- readRDS('AD_GSE147424/pbmc.female.RDS')
MS <- readRDS('MS_GSE193770/pbmc.female.RDS')
pSS <- readRDS('pSS_GSE157278/pbmc.female.RDS')
SLE <- readRDS('SLE_SDY997/pbmc.female.RDS')
UC1 <- readRDS('UC_GSE125527/pbmc.female.RDS')
UC2 <- readRDS('UC_GSE182270/pbmc.female.RDS')

celltype <- commandArgs(trailingOnly=TRUE)
print(celltype)

# Subset each seurat object to the cell type of interest
AD <- subset(AD, cellTypist %in% celltype[[1]], features=rownames(chrX))
MS <- subset(MS, cellTypist %in% celltype[[1]], features=rownames(chrX))
pSS <- subset(pSS, cellTypist %in% celltype[[1]], features=rownames(chrX))
SLE <- subset(SLE, cellTypist %in% celltype[[1]], features=rownames(chrX))
UC1 <- subset(UC1, cellTypist %in% celltype[[1]], features=rownames(chrX))
UC2 <- subset(UC2, cellTypist %in% celltype[[1]], features=rownames(chrX))

# Find overlap between features in each seurat object
feature.list <- list(AD=rownames(AD), MS=rownames(MS), pSS=rownames(pSS), SLE=rownames(SLE), UC1=rownames(UC1), UC2=rownames(UC2))
common.features  <- names(which(table(unlist(feature.list)) == 6))
write.table(common.features, paste0('../common.features/', gsub(' |-|/|', '.', celltype), '.chrX.txt'), sep='\t', quote=F, row.names=F, col.names=F)
