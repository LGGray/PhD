# Script for building Seurat object from SDY998 data

library(Seurat)
library(ddqcR)
library(SeuratDisk)

setwd('datasets/SDY998')

exprMat <- read.delim('celseq_matrix_ru1_reads.tsv', header=T, row.names = 1)
exprMat[is.na(exprMat)] <- 0
meta <- read.delim('celseq_meta.tsv')
pbmc <- CreateSeuratObject(counts = exprMat)
pbmc@meta.data <- cbind(pbmc@meta.data, meta)

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
saveRDS(pbmc, 'pbmc.RDS')

pdf('DimPlot.pdf')
DimPlot(pbmc, label = TRUE, reduction='ref.umap', repel=T) + NoLegend()
dev.off()



