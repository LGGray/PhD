library(Seurat)

# Read in Seurat object
pbmc <- readRDS('pbmc.RDS')
# Read in magic matrix
magic <- read.csv('exp_magic.csv', row.names=1)
# Convert to Seurat assay
magic_assay <- CreateAssayObject(counts=magic)
# Add assay to Seurat object
pbmc[['magic']] <- magic_assay

exp <- GetAssayData(pbmc, slot='counts')
cor_matrix <- cor(magic[rownames(magic) %in% high_cor,], method='spearman')
# identify pairs of features with correlation coefficient > 0.7
high_cor <- which(abs(cor_matrix) > 0.98 & upper.tri(cor_matrix), arr.ind = TRUE)

high_cor <- rownames(exp)[high_cor[,1]]

