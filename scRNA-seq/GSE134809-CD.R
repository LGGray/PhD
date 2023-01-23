library(GEOquery)
library(Seurat)
library(SeuratDisk)
library(ddqcR)

#options(download.file.method.GEOquery = "wget")

#eList2 <- getGEOSuppFiles('GSE134809')
#untar('GSE134809_RAW.tar')

load('~/datasets/XCI/chrY.Rdata')

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/CD_GSE134809')

samples <- '_67|_126|_127|_134|_157|_179|_185|_191|_194'
genes.path <- list.files('.', pattern = 'genes.tsv.gz', recursive = T)
genes.path <- genes.path[grep(samples, genes.path)]
matrix.path <- list.files('.', pattern = 'matrix.mtx.gz', recursive = T)
matrix.path <- matrix.path[grep(samples, matrix.path)]
barcodes.path <- list.files('.', pattern='barcodes.tsv.gz', recursive = T)
barcodes.path <- barcodes.path[grep(samples, barcodes.path)]

# seurat.list <- list()
# seurat.list <- for(i in 1:length(genes.path)){
#   exp.mtx <- ReadMtx(mtx = matrix.path[i], 
#                      features=genes.path[i], 
#                      cells=barcodes.path[i])
#   seurat.list[[i]] <- CreateSeuratObject(counts=exp.mtx)
# }
# Reading in files seperately
exp.mtx <- ReadMtx(mtx = matrix.path[1], 
                      features=genes.path[1], 
                      cells=barcodes.path[1])
GSM4761136 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[2], 
                   features=genes.path[2], 
                   cells=barcodes.path[2])
GSM4761137 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[3], 
                   features=genes.path[3], 
                   cells=barcodes.path[3])
GSM4761138 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[4], 
                   features=genes.path[4], 
                   cells=barcodes.path[4])
GSM4761139 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[5], 
                   features=genes.path[5], 
                   cells=barcodes.path[5])
GSM4761140 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[6], 
                   features=genes.path[6], 
                   cells=barcodes.path[6])
GSM4761141 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[7], 
                   features=genes.path[7], 
                   cells=barcodes.path[7])
GSM4761142 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[8], 
                   features=genes.path[8], 
                   cells=barcodes.path[8])
GSM4761143 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[9], 
                   features=genes.path[9], 
                   cells=barcodes.path[9])
GSM4761144 <- CreateSeuratObject(counts=exp.mtx)

# Merge files together
pbmc <- merge(GSM4761136, y = c(GSM4761137, GSM4761138, GSM4761139,
                                GSM4761140, GSM4761141, GSM4761142,
                                GSM4761143, GSM4761144), 
              add.cell.ids = c('GSM4761136', 'GSM4761137', 'GSM4761138', 
                               'GSM4761139','GSM4761140', 'GSM4761141', 
                               'GSM4761142','GSM4761143', 'GSM4761144'), project = "CD")
# Add individual  and condition to metadata
pbmc$individual <- gsub('_.+', '', rownames(pbmc@meta.data))
pbmc$condition <- 'CD'

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
female <- unique(which.min(table(km.res$cluster)))
sex.list <- ifelse(km.res$cluster == female, 'F', 'M')
# Correcting annotation
sex.list[1] <- 'M'

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