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

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/UC_GSE182270')

files <- list.files(recursive = T, pattern='.tsv|mtx')
files <- files[grep('Raw', files)]

individuals <- dirname(files) %>% 
  gsub('_Grch38/Raw', '', .) %>%
  unique()

expr.list <- lapply(individuals, function(x){
  tmp <- unique(dirname(files[grep(x, files)]))
  Read10X(tmp)
})

pbmc.list <- list()
for(i in 1:length(expr.list)){
  pbmc.list[i] <- CreateSeuratObject(counts = expr.list[[i]])
}
pbmc <- merge(pbmc.list[[1]], c(pbmc.list[[2]], pbmc.list[[3]], pbmc.list[[4]],
                              pbmc.list[[5]], pbmc.list[[6]], pbmc.list[[7]],
                              pbmc.list[[8]]),
              add.cell.ids = individuals)
pbmc$individual <- gsub('_[ACTG]+-.', '', colnames(pbmc))
pbmc$condition <- gsub('\\d.+', '', pbmc$individual)
pbmc$condition <- ifelse(pbmc$condition == 'UC', 'disease', 'control')

# Read in features, barcodes and matrix in working directory
pbmc.data <- Read10X('.')
# Creat Seurat object
pbmc <- CreateSeuratObject(counts=pbmc.data)

# Add Sample information to object
metadata <- read.delim('cell_batch.tsv.gz')
metadata <- metadata[,-3]
pbmc@meta.data <- cbind(pbmc@meta.data, metadata)

# Remove obvious bad quality cells
pbmc.pre <- pbmc
pbmc <- initialQC(pbmc)
# Return dataframe of filtering statistics
pdf('ddqc.plot.pdf')
df.qc <- ddqc.metrics(pbmc, )
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

# Cautious # PC = 18
pbmc <- FindNeighbors(pbmc,dims=1:18)
# Number of clusters in original paper = 47
library("leiden")
leiden_clustering <- leiden(pbmc@graphs$RNA_snn, resolution_parameter = 2)
pbmc@meta.data$leiden_clustering <- leiden_clustering
Idents(pbmc) <- 'leiden_clustering'
pbmc <- RunUMAP(pbmc, dims = 1:18)

pdf('seurat.clusters.DimPlot.pdf')
DimPlot(pbmc, reduction='umap', label=F)
dev.off()

# Remove ambient RNA with decontX
pbmc.raw <- CreateSeuratObject(counts=pbmc.data)
pbmc.raw <- GetAssayData(pbmc.raw, slot = 'counts')
pbmc.expr <- GetAssayData(pbmc, slot = 'counts')
decontaminate <- decontX(pbmc.expr, background=pbmc.raw, z=pbmc$leiden_clustering)
pbmc[["decontXcounts"]] <- CreateAssayObject(counts = decontaminate$decontXcounts)
DefaultAssay(pbmc) <- "decontXcounts"

# Normalise decontX data with Delta method-based variance stabilizing
exp.matrix <- GetAssayData(pbmc, assay = 'decontXcounts', slot = 'counts')
exp.matrix.transformed <- acosh_transform(exp.matrix)

# Add transformed data to Seurat object
pbmc <- SetAssayData(object=pbmc, assay='decontXcounts', slot = 'data', new.data=exp.matrix.transformed)

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
DimPlot(pbmc, reduction='umap', label=TRUE)
dev.off()

# Save complete object
saveRDS(pbmc, 'pbmc.RDS')

# Predicting sex
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/predict.sex.R')

pbmc <- predict.sex(pbmc)

# Subset for male and female
pbmc.male <- subset(pbmc, sex == 'M')
saveRDS(pbmc.male, 'pbmc.male.RDS')

pbmc.female <- subset(pbmc, sex == 'F')
saveRDS(pbmc.female, 'pbmc.female.RDS')

pdf('seurat.cellTypist.pdf', width=15, height=15)
DimPlot(pbmc.female, reduction='umap', label=TRUE)
dev.off()