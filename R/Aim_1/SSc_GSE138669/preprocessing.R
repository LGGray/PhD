library(Seurat)
library(magrittr)
library(ddqcR)
library(reticulate)
library(celda)
library(BiocParallel)
library(scDblFinder)
library(transformGamPoi)
library(iasva)
library(sva)
library(irlba)
library(SummarizedExperiment)

# options(future.globals.maxSize = 891289600)

# scipy_sparse <- import("scipy.sparse")
# mtx <- scipy_sparse$load_npz("exp_counts.npz")
# features <- read.delim('features.tsv.gz', header=F)
# barcodes <- read.delim('barcodes.tsv.gz', header=F)
# colnames(mtx) <- features$V2
# rownames(mtx) <- barcodes$V1
# mtx <- t(mtx)
# Matrix::writeMM(mtx, 'matrix.mtx')

# Read in features, barcodes and matrix in wd
pbmc.data <- Read10X('.')
# Creat Seurat object
pbmc <- CreateSeuratObject(counts=pbmc.data)

# # Create Seurat object
# mtx <- Matrix::readMM('matrix.mtx.gz')
# rownames(mtx) <- features$V2
# colnames(mtx) <- barcodes$V1
# pbmc <- CreateSeuratObject(counts = mtx)

metadata <- read.delim('cell_batch.tsv.gz', row.names = 1)
condition <- read.delim('metadata.txt', sep=' ')
condition <- merge(metadata, condition, by='individual')
metadata$condition <- condition$condition
# Add metadata to Seurat object
pbmc@meta.data <- cbind(pbmc@meta.data, metadata[,c(1,3)])

# Remove obvious bad quality cells
pbmc <- initialQC(pbmc)
# Return dataframe of filtering statistics
pdf('ddqc.plot.pdf')
df.qc <- ddqc.metrics(pbmc)
dev.off()
# Filter out the cells
pbmc <- filterData(pbmc, df.qc)

# Remove ambient RNA with decontX
decontaminate <- decontX(GetAssayData(pbmc, slot = 'counts'))
pbmc[["decontXcounts"]] <- CreateAssayObject(counts = decontaminate$decontXcounts)
DefaultAssay(pbmc) <- "decontXcounts"

# # filter out cells with sum less than 500 for scDblFinder
# pbmc <- subset(pbmc, subset = nCount_RNA > 500)

# remove doublets
sce <- as.SingleCellExperiment(pbmc)
sce$cluster <- fastcluster(sce)
sce <- scDblFinder(sce, samples="individual", clusters=TRUE, BPPARAM=MulticoreParam(3))
pbmc$scDblFinder <- sce$scDblFinder.class
pbmc <- subset(pbmc, scDblFinder == 'singlet')

# Normalise data with Delta method-based variance stabilizing
exp.matrix <- GetAssayData(pbmc, slot = 'counts')
exp.matrix.transformed <- acosh_transform(exp.matrix)

# Add transformed data to Seurat object
pbmc <- SetAssayData(object=pbmc, slot = 'counts', new.data=exp.matrix.transformed)

# Cell type clustering and inspection of markers
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)

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

pbmc <- FindNeighbors(pbmc, dims=1:pcs)
pbmc <- FindClusters(pbmc, resolution=1)
pbmc <- RunUMAP(pbmc, dims = 1:pcs)

pdf('seurat.clusters.DimPlot.pdf')
DimPlot(pbmc, reduction='umap')
dev.off()

# Subset for highly variable genes
HVG <- subset(pbmc, features = VariableFeatures(pbmc))
# Calculate geometric library size
geo_lib_size <- colSums(log(HVG@assays$RNA@data +1))
# Run IA-SVA
set.seed(100)
individual <- pbmc$individual
mod <- model.matrix(~individual + geo_lib_size)
# create a SummarizedExperiment class
sce <- SummarizedExperiment(assay=as.matrix(HVG@assays$RNA@data))
iasva.res <- iasva(sce, mod[, -1], num.sv = 2)
saveRDS(iasva.res, 'iasva.res.RDS')

# Save matrix file for downstream cellTypist analysis
mtx <- as.matrix(GetAssayData(pbmc, slot = 'data'))
write.csv(mtx, 'raw.counts.csv')

# Save unlabelled Seurat object
saveRDS(pbmc, 'pbmc.unlabelled.RDS')
