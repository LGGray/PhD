# Script for building and processing Seurat object from UC_GSE182270 data

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

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/UC_GSE182270')

# files <- list.files(recursive = T, pattern='.tsv|mtx')
# files <- files[grep('Raw', files)]

# individuals <- dirname(files) %>% 
#   gsub('_Grch38/Raw', '', .) %>%
#   unique()

# expr.list <- lapply(individuals, function(x){
#   tmp <- unique(dirname(files[grep(x, files)]))
#   Read10X(tmp)
# })

# pbmc.list <- list()
# for(i in 1:length(expr.list)){
#   pbmc.list[i] <- CreateSeuratObject(counts = expr.list[[i]])
# }
# pbmc <- merge(pbmc.list[[1]], c(pbmc.list[[2]], pbmc.list[[3]], pbmc.list[[4]],
#                               pbmc.list[[5]], pbmc.list[[6]], pbmc.list[[7]],
#                               pbmc.list[[8]]),
#               add.cell.ids = individuals)
# pbmc$individual <- gsub('_[ACTG]+-.', '', colnames(pbmc))
# pbmc$condition <- gsub('\\d.+', '', pbmc$individual)
# pbmc$condition <- ifelse(pbmc$condition == 'UC', 'disease', 'control')



# # Export matrix for DoubletDetector
# library(Matrix)
# mtx <- pbmc@assays$RNA@data
# writeMM(mtx, 'matrix.mtx')
# features <- data.frame(rep(NA, nrow(mtx)), rownames(mtx), rep('Gene Expression', nrow(mtx)))
# write.table(features, 'features.tsv', row.names=F, col.names=F, quote=F, sep='\t')
# barcodes <- colnames(mtx)
# write.table(barcodes, 'barcodes.tsv', row.names=F, col.names=F, quote=F, sep='\t')
# metadata <- pbmc@meta.data[,-c(1:3)]
# metadata$batch <- metadata$individual
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

