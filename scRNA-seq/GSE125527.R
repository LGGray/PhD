library(Seurat)
library(ddqcR)
library(SeuratDisk)

load('~/datasets/XCI/chrY.Rdata')

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/UC_GSE125527')

counts <- read.csv("GSE125527_UMI_cell_table_sparse.csv", header=T)
features <- read.csv("GSE125527_gene_id_rownames.csv", header=F)
cell.names <- read.csv("GSE125527_cell_id_colnames.csv", header=F)
metadata <- read.csv("GSE125527_cell_metadata.csv", row.names = 1, header=T)
colnames(metadata)[c(2,4)] <- c('individual', 'condition') 
metadata$condition <- ifelse(metadata$condition == 'healthy', 'HC', 'UC')

counts.matrix <- matrix(nrow=nrow(features), ncol=nrow(cell.names))
for (i in 1:nrow(counts)){
  counts.matrix[counts[i,1], counts[i,2]] <- counts[i,3]
}
# Replace NA with 0
counts.matrix[is.na(counts.matrix)] <- 0

rownames(counts.matrix) <- features$V1
colnames(counts.matrix) <- cell.names$V1

pbmc <- CreateSeuratObject(counts.matrix)
pbmc@meta.data <- cbind(pbmc@meta.data, metadata)

# Select PBMC data
pbmc <- subset(pbmc, tissue_assignment == 'PBMC')

# Remove obvious bad quality cells
pbmc <- initialQC(pbmc)

# Return dataframe of filtering statistics
pdf('ddqc.plot.pdf')
df.qc <- ddqc.metrics(pbmc)
dev.off()

# Filter out the cells
pbmc <- filterData(pbmc, df.qc)

# SCTransform
pbmc <- SCTransform(pbmc, verbose = FALSE)

# Read in PBMC reference dataset
reference <- LoadH5Seurat("~/azimuth.reference/pbmc_multimodal.h5seurat")
# Find anchors between reference and query
anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)
# Transfer cell type labels and protein data from reference to query
pbmc <- MapQuery(
  anchorset = anchors,
  query = pbmc,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

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
female <- unique(which.max(table(km.res$cluster)))
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