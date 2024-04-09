# Script for building Seurat object from GSE193770 data
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

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/MS_GSE193770')

# files <- list.files(pattern='counts.txt', recursive = T, full.names = T)
# data <- lapply(files, function(x) read.delim(x))
# names(data) <- lapply(files, function(x){
#   unlist(strsplit(basename(x), '_'))[2]
# })

# # Create metadata
# metadata <- lapply(seq_along(data), function(x){
#   data.frame(cell_id=colnames(data[[x]]),
#              batch=rep(names(data)[[x]], ncol(data[[x]])),
#              condition=rep(gsub('\\d', '', names(data)[[x]]), ncol(data[[x]])))
# })

# metadata <- do.call("rbind", metadata)
# metadata$condition <- ifelse(metadata$condition == 'MS', 'disease', 'control')
# write.table(metadata, 'cell_batch.tsv', row.names=F, col.names=T, quote=F, sep='\t')

# # Bind columns across samples
# exprMat <- do.call("cbind", data)
# # Edit cell names
# colnames(exprMat) <- gsub('MS[0-9]+\\.|HC[0-9]+\\.', '', colnames(exprMat))

# # Save matrix for DoubletDetector
# library(Matrix)
# mtx <- Matrix(as.matrix(exprMat), sparse=T)
# writeMM(mtx, 'matrix.mtx')
# features <- data.frame(rep(NA, nrow(exprMat)), rownames(exprMat), rep('Gene Expression', nrow(exprMat)))
# write.table(features, 'features.tsv', row.names=F, col.names=F, quote=F, sep='\t')
# barcodes <- colnames(exprMat)
# write.table(barcodes, 'barcodes.tsv', row.names=F, col.names=F, quote=F, sep='\t')

# Read in features, barcodes and matrix in working directory
pbmc.data <- Read10X('.')
# Creat Seurat object
pbmc.raw <- CreateSeuratObject(counts=pbmc.data)

# Add Sample information to object
metadata <- read.delim('cell_batch.tsv.gz')
colnames(metadata)[2] <- 'individual'
pbmc@meta.data <- cbind(pbmc@meta.data, metadata)

# Check that cell names are equal
all.equal(gsub('\\.[0-9]', '', colnames(pbmc)), gsub('\\.[0-9]', '', pbmc@meta.data$cell_id))

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
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
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

# Find neighbours
pbmc <- FindNeighbors(pbmc, dims=1:11)

# 8 clusters from the paper
library("leiden")
leiden_clustering <- leiden(pbmc@graphs$RNA_snn, resolution_parameter = 0.2)
pbmc@meta.data$leiden_clustering <- leiden_clustering
Idents(pbmc) <- 'leiden_clustering'
pbmc <- RunUMAP(pbmc, dims = 1:11)

pdf('seurat.clusters.DimPlot.pdf')
DimPlot(pbmc, reduction='umap', label=TRUE, raster=FALSE) + NoLegend()
dev.off()

saveRDS(pbmc, 'pbmc.RDS')

# Remove ambient RNA with decontX
pbmc.raw <- GetAssayData(pbmc.raw, slot = 'counts')
pbmc.expr <- GetAssayData(pbmc, slot = 'counts')
decontaminate <- decontX(pbmc.expr, background=pbmc.raw, z=pbmc$leiden_clustering)
pbmc[["decontXcounts"]] <- CreateAssayObject(counts = decontaminate$decontXcounts)
# DefaultAssay(pbmc) <- "decontXcounts"

# Normalise decontX data with Delta method-based variance stabilizing
exp.matrix <- GetAssayData(pbmc, assay = 'decontXcounts', slot = 'counts')
exp.matrix.transformed <- acosh_transform(exp.matrix)

# Add transformed data to Seurat object
pbmc <- SetAssayData(object=pbmc, assay='decontXcounts', slot = 'data', new.data=exp.matrix.transformed)
DefaultAssay(pbmc) <- "decontXcounts"

# Save matrix file for downstream cellTypist analysis
library(data.table)
mtx <- data.frame(GetAssayData(pbmc, assay='decontXcounts', slot = 'counts'))
data.table::fwrite(mtx, file='decontXcounts.counts.csv', row.names=TRUE)

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

# immune_cells <- levels(pbmc)[!(levels(pbmc) %in% c("Megakaryocyte precursor", "Late erythroid", "Megakaryocytes/platelets"))]
# # Remove non-immune cells
# pbmc <- subset(pbmc, cellTypist %in% immune_cells)

### Predict Sex ###
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/predict.sex.R')
pbmc <- predict.sex(pbmc, assay='decontXcounts', slot='data', individual='individual')

# # save object
saveRDS(pbmc, 'pbmc.RDS')

pbmc.female <- subset(pbmc, sex == 'F')
saveRDS(pbmc.female, 'pbmc.female.RDS')