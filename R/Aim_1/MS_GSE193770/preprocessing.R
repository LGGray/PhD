# Script for building Seurat object from GSE193770 data
library(Seurat)
library(magrittr)
library(ddqcR)
library(transformGamPoi)
library(sva)

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

# Read in features, barcodes and matrix in wd
pbmc.data <- Read10X('.')
# Creat Seurat object
pbmc <- CreateSeuratObject(counts=pbmc.data)

# Add Sample information to object
metadata <- read.delim('cell_batch.tsv.gz')
colnames(metadata)[2] <- 'individual'
pbmc@meta.data <- cbind(pbmc@meta.data, metadata)

# Check that cell names are equal
all.equal(gsub('\\.[0-9]', '', colnames(pbmc)), gsub('\\.[0-9]', '', pbmc@meta.data$cell_id))

# Read in DoubletDetector results
dd <- read.delim('DoubletDetection_doublets_singlets.tsv')
dd$CellID <- colnames(pbmc)
dd <- dd[dd$DoubletDetection_DropletType == 'singlet',]
# Remove doublets
pbmc <- subset(pbmc, cells = dd[,1])

# Remove obvious bad quality cells
pbmc <- initialQC(pbmc)

# Return dataframe of filtering statistics
pdf('ddqc.plot.pdf')
df.qc <- ddqc.metrics(pbmc)
dev.off()

# Filter out the cells
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

# Identify batch effects with SVA
# Remove lowly expressed genes
pbmc.HVG <- subset(pbmc, features=VariableFeatures(pbmc))
exp <- pbmc.HVG@assays$RNA@counts
# Full model matrix with variable of interest
mod <- model.matrix(~condition, data=pbmc.HVG@meta.data)
# Null model matrix (include only intercept)
mod0 <- model.matrix(~1, data=pbmc.HVG@meta.data)
# Estimate number of latent factors
n.sv <- num.sv(as.matrix(exp), mod, method='leek')
# Estimate surrogate variables
svseq <- svaseq(as.matrix(exp), mod, mod0, n.sv=n.sv)
save(svseq, file='svaseq.RData')

# Save matrix file for downstream cellTypist analysis
mtx <- as.matrix(GetAssayData(pbmc, slot = 'data'))
write.csv(mtx, 'raw.counts.csv')