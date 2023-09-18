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

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/pSS_GSE157278')

# Read in features, barcodes and matrix in wd
pbmc.data <- Read10X('.')
# Creat Seurat object
pbmc <- CreateSeuratObject(counts=pbmc.data)
# Add Sample information to object
metadata <- read.delim('cell_batch.tsv.gz')
pbmc$individual <- metadata$batch
pbmc$condition <- gsub('-[0-9]', '', metadata$batch)
pbmc$condition <- ifelse(pbmc$condition == 'pSS', 'disease', 'control')

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

# # Harmony integration
# pbmc <- RunHarmony(pbmc, group.by.vars = c('individual'), dims.use = 1:25)

# Being cautious and using 25 PCs
pbmc <- FindNeighbors(pbmc,dims=1:25)
# Number of clusters in original paper = 19
library("leiden")
leiden_clustering <- leiden(pbmc@graphs$RNA_snn, resolution_parameter = 0.3)
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

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
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

# Plot markers before and after decontX
markers <- list(Naive_CD4_T_cell = c('CCR7', 'LTB', 'SELL'),
    Effector_memory_CD4_T_cell = c('CXCR4', 'TNFRSF4', 'CCR6'),
    TH17_cell = c('RORA', 'IL6ST', 'IL17RA'),
    CD4_CTL = c('GZMH', 'GZMA', 'GZMB'),
    CD4_TRAV13_2 = c('TRAV13-2', 'TRBV7-9'),
    Naive_CD8_T_cell = c('CCR7', 'LEF1', 'CD27'),
    CD8_CTL = c('GZMH', 'GZMB', 'ZNF683'),
    MAIT = c('SLC4A10', 'KLRB1', 'ZBTB16'),
    Naive_B_cell = c('CXCR4', 'CD83', 'IGHD'),
    Memory_B_cell = c('CD27'),
    Plasma_cells = c('IGHA1', 'IGHG1', 'IGLC2'),
    pDC = c('CLEC4C', 'GZMB', 'PTPRS', 'IL3RA')
)
pdf('Hong.markers.pdf')
plotDecontXMarkerPercentage(decontaminate$decontXcounts, markers=markers, z=pbmc$leiden_clustering)
dev.off()

bulk.markers <- AggregateExpression(pbmc, features = unique(unlist(markers)), group.by = 'cellTypist', slot = 'counts')$decontXcounts
pdf('Hong.markers.heatmap.pdf')
Heatmap(scale(bulk.markers), col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
column_names_rot = 45, column_names_side = "top", column_dend_side = "bottom")
dev.off()

# Normalise decontX data with Delta method-based variance stabilizing
exp.matrix <- GetAssayData(pbmc, assay = 'decontXcounts', slot = 'counts')
exp.matrix.transformed <- acosh_transform(exp.matrix)

# Add transformed data to Seurat object
pbmc <- SetAssayData(object=pbmc, assay='decontXcounts', slot = 'data', new.data=exp.matrix.transformed)

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

immune_cells <- levels(pbmc)[!(levels(pbmc) %in% c("Megakaryocyte precursor", "Late erythroid", "Megakaryocytes/platelets"))]
# Remove non-immune cells
pbmc <- subset(pbmc, cellTypist %in% immune_cells)

# save objec
saveRDS(pbmc, 'pbmc.female.RDS')