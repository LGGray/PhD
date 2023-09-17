# Script for building and processing Seurat object from UC_GSE125527 data

# Prevent warnings from printing
options(warn=-1)
# Load Libraries
library(Seurat)
library(magrittr)
library(ddqcR)
library(reticulate)
library(celda)
library(BiocParallel)
library(scDblFinder)
library(transformGamPoi)
library(irlba)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/UC_GSE125527')

# counts <- read.csv("GSE125527_UMI_cell_table_sparse.csv", header=T)
# features <- read.csv("GSE125527_gene_id_rownames.csv", header=F)
# cell.names <- read.csv("GSE125527_cell_id_colnames.csv", header=F)
# metadata <- read.csv("GSE125527_cell_metadata.csv", row.names = 1, header=T)
# colnames(metadata)[c(2,4)] <- c('individual', 'condition')
# metadata$batch <- metadata$individual
# metadata$condition <- ifelse(metadata$condition == 'healthy', 'control', 'disease')

# # Create count matrix
# counts.sparse <- sparseMatrix(i=counts$row, j=counts$col, x=counts$value)
# counts.matrix <- as.matrix(counts.sparse)
# rownames(counts.matrix) <- features$V1
# colnames(counts.matrix) <- cell.names$V1

# # Export matrix for DoubletDetector
# library(Matrix)
# mtx <- Matrix(as.matrix(counts.matrix), sparse=T)
# writeMM(mtx, 'matrix.mtx')
# features <- data.frame(rep(NA, nrow(counts.matrix)), rownames(counts.matrix), rep('Gene Expression', nrow(counts.matrix)))
# write.table(features, 'features.tsv', row.names=F, col.names=F, quote=F, sep='\t')
# barcodes <- colnames(counts.matrix)
# write.table(barcodes, 'barcodes.tsv', row.names=F, col.names=F, quote=F, sep='\t')
# write.table(metadata, 'cell_batch.tsv', row.names=F, col.names=T, quote=F, sep='\t')
# # Gzip files
# system('gzip -f matrix.mtx')
# system('gzip -f features.tsv')
# system('gzip -f barcodes.tsv')
# system('gzip -f cell_batch.tsv')

# Read in features, barcodes and matrix in working directory
pbmc.data <- Read10X('.')
# Creat Seurat object
pbmc <- CreateSeuratObject(counts=pbmc.data)

# Add Sample information to object
metadata <- read.delim('cell_batch.tsv.gz')
metadata <- metadata[,-6]
pbmc@meta.data <- cbind(pbmc@meta.data, metadata)

# Remove obvious bad quality cells
pbmc <- initialQC(pbmc)
# Return dataframe of filtering statistics
pdf('ddqc.plot.pdf')
df.qc <- ddqc.metrics(pbmc)
dev.off()
# Filter out the cells
pbmc <- filterData(pbmc, df.qc)

# remove doublets
sce <- as.SingleCellExperiment(pbmc)
sce <- scDblFinder(sce, samples="individual", clusters=TRUE, BPPARAM=MulticoreParam(3))
pbmc$scDblFinder <- sce$scDblFinder.class
pbmc <- subset(pbmc, scDblFinder == 'singlet')

# Normalise data with Delta method-based variance stabilizing
exp.matrix <- GetAssayData(pbmc, slot = 'counts', assay='RNA')
exp.matrix.transformed <- acosh_transform(exp.matrix)

# Add transformed data to Seurat object
pbmc <- SetAssayData(object=pbmc, assay='RNA', slot = 'data', new.data=exp.matrix.transformed)

# Cell type clustering and inspection of markers
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)

# Find the number of PCs to use
# Plot the elbow plot
pdf('seurat.pca.pdf')
ElbowPlot(pbmc, ndims = 50)
dev.off()

# Plot PCA to show individual 
pdf('seurat.pca.individual.pdf')
DimPlot(pbmc, reduction='pca', group.by='individual')
dev.off()
# Or condition effects
pdf('seurat.pca.condition.pdf')
DimPlot(pbmc, reduction='pca', group.by='condition')
dev.off()

pdf('seurat.vlnplot.condition.pdf')
VlnPlot(pbmc, features="PC_1", group.by="condition")
dev.off()

# Determine percent of variation associated with each PC
pct <- pbmc[["pca"]]@stdev / sum(pbmc[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# Minimum of the two calculation
pcs <- min(co1, co2)
print(paste('Selected # PCs', pcs))

# Cautious # PC = 17
# Being cautious and using 25 PCs
pbmc <- FindNeighbors(pbmc,dims=1:25)
# Number of clusters in original paper = 19
library("leiden")
leiden_clustering <- leiden(pbmc@graphs$RNA_snn, resolution_parameter = 2)
pbmc@meta.data$leiden_clustering <- leiden_clustering
Idents(pbmc) <- 'leiden_clustering'
pbmc <- RunUMAP(pbmc, dims = 1:25)

pdf('seurat.clusters.DimPlot.pdf')
DimPlot(pbmc, reduction='umap', label=TRUE)
dev.off()

# Remove ambient RNA with decontX
pbmc.raw <- CreateSeuratObject(counts=pbmc.data)
pbmc.raw <- GetAssayData(pbmc.raw, slot = 'counts')
pbmc.expr <- GetAssayData(pbmc, slot = 'counts')
decontaminate <- decontX(pbmc.expr, background=pbmc.raw, z=pbmc$leiden_clustering)
pbmc[["decontXcounts"]] <- CreateAssayObject(counts = decontaminate$decontXcounts)
DefaultAssay(pbmc) <- "decontXcounts"

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)
pbmc <- ScaleData(pbmc, features=NULL, assay='decontXcounts')

# Plot markers before and after decontX
markers <- list(Tcell_Markers = c("CD3E", "CD3D", "TRAC", "LTB", "SELL", "CCR7"),
    Bcell_Markers = c("CD79A", "CD79B", "MS4A1"),
    Monocyte_Markers = c("S100A8", "S100A9", "LYZ"),
    NKcell_Markers = "GNLY")
pdf('Before.decontX.markers.pdf')
plotDecontXMarkerPercentage(pbmc.expr, markers=markers, z=pbmc$leiden_clustering)
dev.off()

pdf('After.decontX.markers.pdf')
plotDecontXMarkerPercentage(decontaminate$decontXcounts, markers=markers, z=pbmc$leiden_clustering)
dev.off()

# Save matrix file for downstream cellTypist analysis
mtx <- as.matrix(GetAssayData(pbmc, assay='decontXcounts', slot = 'counts'))
write.csv(mtx, 'decontXcounts.counts.csv')

# Save unlabelled Seurat object
saveRDS(pbmc, 'pbmc.unlabelled.RDS')

# Read in cellTypist labels
cellTypist <- read.csv('cellTypist/predicted_labels.csv', row.names=1)
# Add cellTypist labels to Seurat object
pbmc$cellTypist <- cellTypist$majority_voting

# Plot cellTypist labels
Idents(pbmc) <- 'cellTypist'
pdf('seurat.cellTypist.pdf', width=15, height=15)
DimPlot(pbmc, reduction='umap', label=TRUE)
dev.off()

boland.markers <- c('CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD8B', 'CD4','NCAM1', 'TBX21', 
'IFNG', 'GATA3', 'CXCR5', 'BCL6', 'PDCD1', 'FOXP3', 'FCER2', 'ZBTB7B', 'RUNX3',
'CD14', 'ITGAM', 'CEBPA', 'CEBPB', 'SPI1', 'FLT3', 'ITGAX', 'ZBTB46', 'CEACAM8', 
'ZBTB10', 'IL3RA','CLEC4C', 'NRP1', 'NKG7', 'KLRC1', 'KLRC2', 'KLRC3', 'KLRB1', 
'Ly6G6D', 'ADGRE1', 'CD19', 'SDC1', 'MS4A1', 'CD38', 'CD27', 'TNFRSF13C', 'CD84',
'CD86', 'CD24', 'TRAV', 'TRAC', 'TRBV', 'TRBC', 'TRGV', 'TRGC', 'TRDV', 'TRDC', 
'CR2', 'CD22', 'IGHV', 'IGHG', 'IGHA', 'IGHD', 'IGHM', 'IGHE', 'IGLV', 'IGLL',
'IGKV', 'IGLC', 'IGKC')

bulk.markers <- AggregateExpression(pbmc, features = boland.markers, group.by = 'cellTypist', slot = 'counts')$decontXcounts
pdf('Boland.markers.heatmap.pdf')
Heatmap(scale(bulk.markers), col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
column_names_rot = 45, column_names_side = "top", column_dend_side = "bottom")
dev.off()

Idents(pbmc) <- 'cellTypist'

# Save complete object
saveRDS(pbmc, 'pbmc.RDS')

# Predicting sex
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/function')

# Subset for male and female
pbmc.male <- subset(pbmc, sex == 'M')
saveRDS(pbmc.male, 'pbmc.male.RDS')

pbmc.female <- subset(pbmc, sex == 'F')
saveRDS(pbmc.female, 'pbmc.female.RDS')

pdf('seurat.cellTypist.pdf', width=15, height=15)
DimPlot(pbmc.female, reduction='umap', label=TRUE)
dev.off()