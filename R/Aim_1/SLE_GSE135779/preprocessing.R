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

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/SLE_GSE135779/')

files <- list.files('scRNAseq.files/', recursive = T, pattern='.tsv|mtx', full.names=T)

individuals <- strsplit(files[-1], '_') %>% sapply(function(x) x[2]) %>% unique()

meta <- read.csv('samples_data.csv')
features <- read.delim("scRNAseq.files/GSE135779_genes.tsv.gz", header = F, sep = '\t')

seurat_objects <- list()

seurat_objects <- lapply(individuals, function(x){
  tmp <- unique(files[grep(x, files)])
  barcodes <- read.table(tmp[1], header = F, sep = '\t')
  matrix <- readMM(tmp[2])
  colnames(matrix) <- barcodes$V1
  rownames(matrix) <- features$V2

  seurat_object <- CreateSeuratObject(counts = matrix)
  seurat_object$individual <- x
  seurat_object$condition <- meta$condition[meta$individual == x]
  return(seurat_object)
})

pbmc <- seurat_objects[[1]]

# Merge the rest of the Seurat objects
for (i in 2:length(seurat_objects)) {
  pbmc <- merge(pbmc, y = seurat_objects[[i]])
}

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
pbmc <- FindNeighbors(pbmc, dims=1:16)

# 20 clusters from the paper
library("leiden")
leiden_clustering <- leiden(pbmc@graphs$RNA_snn, resolution_parameter = 0.4)
pbmc@meta.data$leiden_clustering <- leiden_clustering
Idents(pbmc) <- 'leiden_clustering'
pbmc <- RunUMAP(pbmc, dims = 1:16)

pdf('seurat.clusters.DimPlot.pdf')
DimPlot(pbmc, reduction='umap', label=TRUE, raster=FALSE) + NoLegend()
dev.off()

saveRDS(pbmc, 'pbmc.RDS')

# Remove ambient RNA with decontX
# pbmc.raw <- CreateSeuratObject(counts=pbmc.raw@assays$RNA@counts)
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
mtx <- as.matrix(GetAssayData(pbmc, assay='decontXcounts', slot = 'counts'))
data.table::fwrite(mtx, file='decontXcounts.counts.csv', row.names=TRUE)

mtx <- data.table::fread('decontXcounts.counts.csv')

# Normalise decontX data with Delta method-based variance stabilizing
exp.matrix <- GetAssayData(pbmc, assay = 'decontXcounts', slot = 'counts')
exp.matrix.transformed <- acosh_transform(exp.matrix)

# Add transformed data to Seurat object
pbmc <- SetAssayData(object=pbmc, assay='decontXcounts', slot = 'data', new.data=exp.matrix.transformed)

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
DimPlot(pbmc, reduction='umap', label=TRUE, raster=F)
dev.off()

# immune_cells <- levels(pbmc)[!(levels(pbmc) %in% c("Megakaryocyte precursor", "Late erythroid", "Megakaryocytes/platelets"))]
# # Remove non-immune cells
# pbmc <- subset(pbmc, cellTypist %in% immune_cells)

### Predict Sex ###
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/predict.sex.R')
pbmc <- predict.sex(pbmc, assay='decontXcounts', slot='data', individual='individual')

pbmc$condition <- gsub('\\d', '', pbmc$condition)

# # save object
# saveRDS(pbmc, 'pbmc.RDS')
SeuratDisk::SaveH5Seurat(pbmc, filename = "pbmc.h5Seurat", overwrite = TRUE)
# pbmc <- LoadH5Seurat("pbmc.h5Seurat")

pbmc.female <- subset(pbmc, sex == 'F')
SeuratDisk::SaveH5Seurat(pbmc.female, filename = "pbmc.female.h5Seurat", overwrite = TRUE)


###### Rename cellTypist labels based on Perez et al. 2022

library(Seurat)
library(dplyr)
library(edgeR)

perez_celltypes <- list.files('../../SLE/pseudobulk/split_1/combined/ensemble/', pattern='.chrX.sav')
perez_celltypes <- gsub('.chrX.sav', '', perez_celltypes)
perez_celltypes <- gsub('_', ' ', perez_celltypes)
perez_celltypes[6] <- "CD4 T cell Effector-Memory"
perez_celltypes[10] <- "CD8 T cell Cytotoxic-GZMH"
perez_celltypes[11] <- "CD8 T cell Cytotoxic-GZMK"
perez_celltypes[20] <- "Non-classical monocyte"

levels(pbmc) %in% perez_celltypes


pbmc <- readRDS('pbmc.female.RDS')

Perez <- pbmc$cellTypist

Perez <- case_when(
  Perez == "NK cells" ~ "natural killer cell",
  Perez == "Naive B cells" ~ "B cell Naive",
  Perez == "B cells" ~ "B cell",
  Perez == "Classical monocytes" ~ "Classical monocyte",
  Perez == "Non-classical monocytes" ~ "Non-classical monocyte",
  Perez == "Tcm/Naive helper T cells" ~ "CD4 T cell Naive",
  Perez == "Tcm/Naive cytotoxic T cells" ~ "CD8 T cell Naive",
  Perez == "Tem/Trm cytotoxic T cells" ~ "CD8 positive, alpha beta T cell",
  Perez == "Tem/Effector helper T cells" ~ "CD4 positive, alpha beta T cell",
  Perez == "Regulatory T cells" ~ "CD4 T cell Treg",
  Perez == "MAIT cells" ~ "CD8 T cell MAIT",
  Perez == "DC1" ~ "Conventional DC",
  Perez == "DC2" ~ "Conventional DC",
  Perez == "CD16+ NK cells" ~ "NK cell Bright",
  Perez == "pDC" ~ "Plasmacytoid DC",
  Perez == "Age-associated B cells" ~ "B cell Atypical",
  Perez == "Plasmablasts" ~ "Plasmablast",
  Perez == "HSC/MPP" ~ "Progenitor cell",
  TRUE ~ Perez
)

pbmc@meta.data$Perez <- Perez

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/Seurat2PB.R')

# Load X chromosome genes
chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name

# Load SLE DisGenNet genes
SLE <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/SLE.tsv')
SLE <- unique(SLE$Gene)

celltypes <- unique(pbmc$Perez)[unique(pbmc$Perez) %in% perez_celltypes]

for(cell in celltypes){
    pbmc.subset <- subset(pbmc, Perez == cell)
    # Pseudobulk by summed expression
    bulk <- AggregateExpression(pbmc.subset, slot='counts', group.by='individual')$RNA
    # Remove lowly expressed genes
    non_zero_genes <- rowSums(bulk) > 0
    bulk <- bulk[non_zero_genes, ]
    # Divide each gene by total number of counts across all genes
    bulk <- apply(bulk, 2, function(x) x / sum(x))

    # Subset to X chromosome genes
    bulk.chrX <- data.frame(t(bulk[rownames(bulk) %in% chrX,]))

    # Create metadata to add to pseudobulk matrix
    meta <- unique(pbmc.subset@meta.data[,c('condition', 'individual')])
    meta$cellCount <- sapply(meta$individual, function(id) sum(pbmc.subset$individual == id))
    meta <- meta[match(rownames(bulk.chrX), meta$individual),]

    # Add metadata to each pseudobulk matrix
    bulk.chrX <- cbind(class=meta$condition, individual=meta$individual, 
    cellCount=meta$cellCount, bulk.chrX)

    # Save pseudobulk matrix
    cell <- gsub("\\+|-| ", "_", cell)
    saveRDS(bulk.chrX, paste0('pseudobulk_update/', cell, '.chrX.RDS'))
}

# Export the pseudobulked expression matrix, subsetted by SLE DisGenNet
for(cell in celltypes){
    pbmc.subset <- subset(pbmc, Perez == cell)
    # Pseudobulk by summed expression
    bulk <- AggregateExpression(pbmc.subset, slot='counts', group.by='individual')$RNA
    # Remove lowly expressed genes
    non_zero_genes <- rowSums(bulk) > 0
    bulk <- bulk[non_zero_genes, ]
    # Divide each gene by total number of counts across all genes
    bulk <- apply(bulk, 2, function(x) x / sum(x))

    # Subset for SLE genes
    bulk.SLE <- data.frame(t(bulk[rownames(bulk) %in% SLE,]))

    # Create metadata to add to pseudobulk matrix
    meta <- unique(pbmc.subset@meta.data[,c('condition', 'individual')])
    meta$cellCount <- sapply(meta$individual, function(id) sum(pbmc.subset$individual == id))
    meta <- meta[match(rownames(bulk.SLE), meta$individual),]

    # Add metadata to each pseudobulk matrix
    bulk.SLE <- cbind(class=meta$condition, individual=meta$individual, 
    cellCount=meta$cellCount, bulk.SLE)

    # Save pseudobulk matrix
    cell <- gsub("\\+|-| ", "_", cell)
    saveRDS(bulk.SLE, paste0('pseudobulk_update/', cell, '.SLE.RDS'))
}

