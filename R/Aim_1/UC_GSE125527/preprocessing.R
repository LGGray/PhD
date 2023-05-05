# Script for building and processing Seurat object from UC_GSE125527 data

# Prevent warnings from printing
options(warn=-1)
# Load Libraries
library(Seurat)
library(magrittr)
library(ddqcR)
library(transformGamPoi)
library(iasva)
library(sva)
library(irlba)
library(SummarizedExperiment)

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
pbmc@meta.data <- cbind(pbmc@meta.data, metadata)

# Check that cell names are equal
all.equal(gsub('\\.[0-9]', '', colnames(pbmc)), gsub('\\.[0-9]', '', pbmc@meta.data$cell_id))

# Read in DoubletDetector results
dd <- read.delim('DoubletDetection_doublets_singlets.tsv')
dd$CellID <- colnames(pbmc)
dd <- dd[dd$DoubletDetection_DropletType == 'singlet',]
# Remove doublets
pbmc <- subset(pbmc, cells = dd[,1])

# Select PBMC data
pbmc <- subset(pbmc, tissue_assignment == 'PBMC')

# Remove obvious bad quality cells
pbmc <- initialQC(pbmc)

# Return dataframe of filtering statistics
pdf('ddqc.plot.pdf')
df.qc <- ddqc.metrics(pbmc)
dev.off()

# Filter out the low quality cells
pbmc <- filterData(pbmc, df.qc)

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
iasva.res <- iasva(sce, mod[, -1], num.sv = 5)
saveRDS(iasva.res, 'iasva.res.RDS')

# Save matrix file for downstream cellTypist analysis
mtx <- as.matrix(GetAssayData(pbmc, slot = 'data'))
write.csv(mtx, 'raw.counts.csv')

# Save unlabelled Seurat object
saveRDS(pbmc, 'pbmc.unlabelled.RDS')