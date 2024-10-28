
# Load required libraries
library(Seurat)
library(NMF)
library(edgeR)
library(ComplexHeatmap)
library(circlize)

if(!dir.exists('NMF') == TRUE) {
  dir.create('NMF')
}

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/Seurat2PB.R')

chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name

pbmc <- readRDS('pbmc.female.RDS')

cell <- levels(pbmc)[commandArgs(trailingOnly = TRUE)[1]]
pbmc.subset <- subset(pbmc, cellTypist == cell)

mtx <- Seurat2PB(pbmc.subset, sample='individual', cluster='condition', assay='RNA')$counts

# remove genes with zero counts
mtx <- mtx[rowSums(mtx) > 0,]

# Set rank and run NMF with 8 cores
rank = 2
nmf_result <- nmf(mtx, rank, nrun=100, method = "brunet", .options = 'P8v')

# Save NMF results
save(nmf_result, file=paste0('NMF/', gsub("/|-| ", "_", cell), '.RData'))

# Extract and save features
s <- extractFeatures(nmf_result)
features <- unique(unlist(lapply(s, function(x) rownames(w)[x])))
write.csv(features, file=paste0('NMF/', gsub("/|-| ", "_", cell), '_features.csv'), row.names=FALSE)

# Subset for chrX genes
features <- features[features %in% chrX]
plot_mtx <- scale(mtx[features,])

pdf(paste0('NMF/', gsub("/|-| ", "_", cell), '_heatmap.pdf'))
Heatmap(plot_mtx, name = "z-score", 
col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), 
show_row_names = TRUE, cluster_columns = FALSE)
dev.off()