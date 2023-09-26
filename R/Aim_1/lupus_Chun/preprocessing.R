library(Seurat)
library(sva)
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


# # scipy_sparse <- import("scipy.sparse")
# # mtx <- scipy_sparse$load_npz("exp_counts.npz")
# # features <- read.delim('features.tsv.gz')
# # barcodes <- read.delim('barcodes.tsv.gz', header=F)
# # colnames(mtx) <- features$feature_name
# # rownames(mtx) <- barcodes$V1
# # mtx <- t(mtx)
# # Matrix::writeMM(mtx, 'matrix.mtx')

# # Read in features, barcodes and matrix in wd
# pbmc.data <- Read10X('.')
# # Creat Seurat object
# pbmc <- CreateSeuratObject(counts=pbmc.data)

# # # Create Seurat object
# # mtx <- Matrix::readMM('matrix.mtx.gz')
# # rownames(mtx) <- features$V2
# # colnames(mtx) <- barcodes$V1
# # pbmc <- CreateSeuratObject(counts = mtx)

# metadata <- read.delim('cell_batch.tsv.gz', row.names = 1)
# colnames(metadata)[26] <- 'condition'
# metadata$condition <- ifelse(metadata$condition == 'systemic lupus erythematosus', 'disease', 'control')
# colnames(metadata)[16] <- 'individual'
# metadata$sex <- ifelse(metadata$sex == 'female', 'F', 'M')
# colnames(metadata)[31] <- 'age'
# metadata$age <- as.numeric(gsub('-year-old human stage', '', metadata$age))
# # Add metadata to Seurat object
# pbmc@meta.data <- cbind(pbmc@meta.data, metadata)

# ### Remove flare and treated samples ###
# pbmc <- subset(pbmc, disease_state %in% c('na', 'managed'))

# # # Adjust counts for batch effects
# # adjusted_counts <- ComBat_seq(exp.matrix, batch=pbmc$Processing_Cohort, group=pbmc$condition)
# # pbmc <- SetAssayData(object=pbmc, assay='RNA', slot = 'counts', new.data=adjusted_counts)

# # # Remove obvious bad quality cells
# # pbmc <- initialQC(pbmc)
# # # Return dataframe of filtering statistics
# # pdf('ddqc.plot.pdf')
# # df.qc <- ddqc.metrics(pbmc)
# # dev.off()
# # # Filter out the cells
# # pbmc <- filterData(pbmc, df.qc)

# # Normalise data with Delta method-based variance stabilizing
# exp.matrix <- GetAssayData(pbmc, slot = 'counts', assay='RNA')
# exp.matrix.transformed <- acosh_transform(exp.matrix)

# # Add transformed data to Seurat object
# pbmc <- SetAssayData(object=pbmc, assay='RNA', slot = 'data', new.data=exp.matrix.transformed)

# # Cell type clustering and inspection of markers
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# pbmc <- ScaleData(pbmc)
# pbmc <- RunPCA(pbmc)

# # Find the number of PCs to use
# # Plot the elbow plot
# pdf('seurat.pca.pdf')
# ElbowPlot(pbmc, ndims = 50)
# dev.off()

# # Determine percent of variation associated with each PC
# pct <- pbmc[["pca"]]@stdev / sum(pbmc[["pca"]]@stdev) * 100
# # Calculate cumulative percents for each PC
# cumu <- cumsum(pct)
# # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
# co1 <- which(cumu > 90 & pct < 5)[1]
# # Determine the difference between variation of PC and subsequent PC
# co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# # Minimum of the two calculation
# pcs <- min(co1, co2)
# print(paste('Selected # PCs', pcs))

pbmc <- readRDS('pbmc.unlabelled.RDS')

# Being cautious and using 20 PCs
pbmc <- FindNeighbors(pbmc,dims=1:20)
# Number of clusters in original paper = 23
library("leiden")
leiden_clustering <- leiden(pbmc@graphs$RNA_snn, resolution_parameter = 3)
pbmc@meta.data$leiden_clustering <- leiden_clustering
Idents(pbmc) <- 'leiden_clustering'
pbmc <- RunUMAP(pbmc, dims = 1:20)

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

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features=NULL, assay='decontXcounts')

# Save matrix file for downstream cellTypist analysis
mtx <- data.frame(GetAssayData(pbmc, assay='decontXcounts', slot = 'counts'))
data.table::fwrite(mtx, 'decontXcounts.counts.csv', row.names=T)

# Save unlabelled Seurat object
saveRDS(pbmc, 'pbmc.unlabelled.RDS')