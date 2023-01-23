library(Seurat)
library(ddqcR)
library(SeuratDisk)

load('~/datasets/XCI/chrY.Rdata')

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/T1D_parsebio')

# Convert("T1D.h5ad", dest="h5seurat", overwrite=F)
pbmc <- LoadH5Seurat("T1D.h5seurat")
metadata <- read.csv('T1D.metadata.csv', row.names=1)
colnames(metadata)[1] <- 'individual'
metadata$condition <- gsub('_[0-9]+', '', metadata$individual)
expr <- GetAssayData(pbmc, assay='RNA')
pbmc <- CreateSeuratObject(expr)
pbmc@meta.data <- cbind(pbmc@meta.data, metadata)

# Remove obvious bad quality cells
pbmc <- initialQC(pbmc)

# Perform the filtering * Based on parsebio website *
pbmc <- subset(pbmc, subset = nFeature_RNA < 5000 & nCount_RNA < 20000 & percent.mt < 15)

# Visualize QC metrics as a violin plot
pdf('QC.Vlnplot.pdf')
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# SCTransform
pbmc <- SCTransform(pbmc, verbose = FALSE)

# Read in PBMC reference dataset
reference <- LoadH5Seurat("~/azimuth.reference/pbmc_multimodal.h5seurat")
# Find anchors between reference and query
anchors <- FindTransferAnchors(reference = reference,
                               query = pbmc,
                               normalization.method = "SCT",
                               reference.reduction = "spca",
                               dims = 1:50)
# Transfer cell type labels and protein data from reference to query
pbmc <- MapQuery(anchorset = anchors,
                 query = pbmc,
                 reference = reference,
                 refdata = list(celltype.l1 = "celltype.l1",
                                celltype.l2 = "celltype.l2",
                                predicted_ADT = "ADT"),
                 reference.reduction = "spca",
                 reduction.model = "wnn.umap")

Idents(pbmc) <- 'predicted.celltype.l2'

pdf('DimPlot.all.pdf')
DimPlot(pbmc, label = TRUE, reduction='ref.umap', repel=T) + NoLegend()
dev.off()

# Determine sex of individuals by expression of chrY and XIST
set.seed(42)
exp <- AverageExpression(pbmc, assays='SCT', features=c('XIST', rownames(chrY)), group.by='individual')
pca <- prcomp(t(exp[[1]]), scale=T)
pdf('sex.PCA.pdf')
plot(pca$x[,1], pca$x[,2])
text(pca$x[,1], pca$x[,2], labels = colnames(exp[[1]]))
dev.off()

# K-means clustering on the PCA data *Males more common in study*
pc <- pca$x
km.res <- kmeans(pc[,1:2], centers = 2, nstart = 50)
female <- unique(which.min(table(km.res$cluster)))
sex.list <- ifelse(km.res$cluster == female, 'F', 'M')

# Add sex to metadata
pbmc$sex <- sex.list[pbmc$individual]

# Output all cells
saveRDS(pbmc, 'pbmc.RDS')
# Subset for females and output 
pbmc.female <- subset(pbmc, sex == 'F')

pdf('DimPlot.female.pdf')
DimPlot(pbmc.female, label = TRUE, reduction='ref.umap', repel=T) + NoLegend()
dev.off()

DefaultAssay(pbmc) <- 'RNA'
saveRDS(pbmc.female, 'pbmc.female.RDS')
