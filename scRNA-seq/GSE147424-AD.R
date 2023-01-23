library(GEOquery)
library(Seurat)
library(SeuratDisk)
library(ddqcR)

load('~/datasets/XCI/chrY.Rdata')

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/AD_GSE147424')

options(download.file.method.GEOquery = "wget")

eList2 <- getGEOSuppFiles('GSE147424')
untar('GSE147424_RAW.tar')

# Read in files
files <- list.files(pattern='.gz')
seurat.list <- lapply(files, function(x) {
  exp <- read.delim(x, sep=',', row.names = 1)
  return(CreateSeuratObject(counts=exp))
})

pbmc <- merge(seurat.list[[1]], y=c(seurat.list[[2]], seurat.list[[3]], seurat.list[[4]], seurat.list[[5]], 
                                    seurat.list[[6]], seurat.list[[7]], seurat.list[[8]], seurat.list[[9]], 
                                    seurat.list[[10]], seurat.list[[11]], seurat.list[[12]], seurat.list[[13]], 
                                    seurat.list[[14]], seurat.list[[15]], seurat.list[[16]], seurat.list[[17]]),
              add.cell.ids = gsub('_.+', '', files))
# Adding metadata
pbmc$individual <- gsub('_.+', '', rownames(pbmc@meta.data))
condition.list <- c(GSM4430459='LS', GSM4430460='LS', GSM4430461='NL', GSM4430462='HC', GSM4430463='LS',
                       GSM4430464='HC', GSM4430465='LS', GSM4430466='HC', GSM4430467='HC', GSM4430468='HC',
                       GSM4430469='NL', GSM4430470='HC', GSM4430471='HC', GSM4430472='NL', GSM4430473='NL',
                       GSM4430474='NL', GSM4430475='HC')
pbmc$condition <- condition.list[pbmc$individual]

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

# Clustering cells
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), ass)
pdf('elbowplot.pdf')
ElbowPlot(pbmc)
dev.off()

pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.5)

all.markers <- FindAllMarkers(pbmc)
all.markers.split <- split(all.markers, all.markers$cluster)
all.markers.split 
save(all.markers.split, file='cluster.markers.RData')



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