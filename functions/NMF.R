
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

cell <- levels(pbmc)[as.numeric(commandArgs(trailingOnly = TRUE)[1])]
print(cell)
pbmc.subset <- subset(pbmc, cellTypist == cell)

mtx <- Seurat2PB(pbmc.subset, sample='individual', cluster='condition', assay='RNA')$counts

# remove genes with zero counts
mtx <- mtx[rowSums(mtx) > 0,]

# Set rank and run NMF with 8 cores
rank = 2
nmf_result <- nmf(mtx, rank, nrun=100, method = "brunet", .options = 'Pv')

# Save NMF results
save(nmf_result, file=paste0('NMF/', gsub("/|-| ", "_", cell), '.RData'))

# Extract and save features
s <- extractFeatures(nmf_result)
features <- unique(unlist(lapply(s, function(x) rownames(w)[x])))
write.csv(features, file=paste0('NMF/', gsub("/|-| ", "_", cell), '_features.csv'), row.names=FALSE)

# Subset for chrX genes
features <- features[features %in% chrX]
plot_mtx <- scale(mtx[rownames(mtx) %in% features,])
colnames(plot_mtx) <- gsub('.+_cluster', '', colnames(plot_mtx))

pdf(paste0('NMF/', gsub("/|-| ", "_", cell), '_heatmap.pdf'))
# ha <- HeatmapAnnotation(df = data.frame(condition = colnames(plot_mtx)), 
# col = list(condition = c(control = "black", disease = "red")))
Heatmap(plot_mtx, name = "z-score", 
col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), 
show_row_names = TRUE, cluster_columns = FALSE, 
show_column_names = FALSE,
column_split = colnames(plot_mtx))
dev.off()