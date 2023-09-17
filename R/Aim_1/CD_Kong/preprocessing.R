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

# args <- commandArgs(trailingOnly = TRUE)

# setwd(paste0('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/CD_Kong/', args[1]))

# scipy_sparse <- import("scipy.sparse")
# mtx <- scipy_sparse$load_npz("exp_counts.npz")
# features <- read.delim('features.tsv.gz')
# barcodes <- read.delim('barcodes.tsv.gz', header=F)
# colnames(mtx) <- features$feature_name
# rownames(mtx) <- barcodes$V1
# mtx <- t(mtx)
# Matrix::writeMM(mtx, 'matrix.mtx')

#gzip matrix.mtx
# R.utils::gzip('matrix.mtx', overwrite=T)

# Read in features, barcodes and matrix in wd
pbmc.data <- Read10X('.')
# Creat Seurat object
pbmc <- CreateSeuratObject(counts=pbmc.data)
# Add Sample information to object
metadata <- read.delim('cell_batch.tsv.gz')
pbmc@meta.data <- cbind(pbmc@meta.data, metadata)
colnames(pbmc@meta.data)[9] <- 'individual'
colnames(pbmc@meta.data)[25] <- 'condition'
pbmc$condition <- ifelse(pbmc$condition == 'Crohn disease', 'disease', 'control')

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
sce <- scDblFinder(sce, samples="individual", clusters='cell_type', BPPARAM=MulticoreParam(3))
pbmc$scDblFinder <- sce$scDblFinder.class
pbmc <- subset(pbmc, scDblFinder == 'singlet')

# Normalise data with Delta method-based variance stabilizing
exp.matrix <- GetAssayData(pbmc, slot = 'counts')
exp.matrix.transformed <- acosh_transform(exp.matrix)

# Add transformed data to Seurat object
pbmc <- SetAssayData(object=pbmc, slot = 'data', new.data=exp.matrix.transformed)

# Cell type clustering and inspection of markers
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)

# Plot PCA to show individual 
pdf('seurat.pca.individual.pdf')
DimPlot(pbmc, reduction='pca', group.by='individual', raster=F)
dev.off()
# Or condition effects
pdf('seurat.pca.condition.pdf')
DimPlot(pbmc, reduction='pca', group.by='condition', raster=F)
dev.off()
# Region effects
pdf('seurat.pca.batch.pdf')
DimPlot(pbmc, reduction='pca', group.by='Type', raster=F) + NoLegend()
dev.off()

pdf('seurat.vlnplot.condition.pdf')
VlnPlot(pbmc, features="PC_1", group.by="Layer")
dev.off()


# Find the number of PCs to use
# Plot the elbow plot
pdf('seurat.pca.pdf')
ElbowPlot(pbmc, ndims = 50)
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

# Cautious # PC = 20
pbmc <- FindNeighbors(pbmc,dims=1:20)
# Number of clusters in original paper = 19
library("leiden")
leiden_clustering <- leiden(pbmc@graphs$RNA_snn, resolution_parameter = 1.3)
pbmc@meta.data$leiden_clustering <- leiden_clustering
Idents(pbmc) <- 'leiden_clustering'
pbmc <- RunUMAP(pbmc, dims = 1:20)

pdf('seurat.clusters.DimPlot.pdf')
DimPlot(pbmc, reduction='umap', raster=FALSE, label=TRUE)
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
DefaultAssay(pbmc) <- "decontXcounts"

# Normalise decontX data with Delta method-based variance stabilizing
exp.matrix <- GetAssayData(pbmc, assay = 'decontXcounts', slot = 'counts')
exp.matrix.transformed <- acosh_transform(exp.matrix)

# Add transformed data to Seurat object
pbmc <- SetAssayData(object=pbmc, assay='decontXcounts', slot = 'data', new.data=exp.matrix.transformed)


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
mtx <- data.frame(GetAssayData(pbmc, assay='decontXcounts', slot = 'counts'))
data.table::fwrite(mtx, 'decontXcounts.counts.csv', row.names=T)

# Save unlabelled Seurat object
saveRDS(pbmc, 'pbmc.unlabelled.RDS')

# Read in cellTypist labels
cellTypist <- read.csv('cellTypist/predicted_labels.csv', row.names=1)
# Add cellTypist labels to Seurat object
pbmc$cellTypist <- cellTypist$majority_voting

# Plot cellTypist labels
Idents(pbmc) <- 'cellTypist'
pdf('seurat.cellTypist.pdf', width=15, height=15)
DimPlot(pbmc, reduction='umap', raster=FALSE, label=TRUE)
dev.off()

# Subset for male and female and save object
pbmc.male <- subset(pbmc, sex == 'male')
saveRDS(pbmc.male, 'pbmc.male.RDS')

pbmc.female <- subset(pbmc, sex =='female')
saveRDS(pbmc.female, 'pbmc.female.RDS')

pdf('seurat.cellTypist.pdf', width=15, height=15)
DimPlot(pbmc.female, reduction='umap', raster=FALSE, label=TRUE)
dev.off()