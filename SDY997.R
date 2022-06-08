# Script for building Seurat object from SDY997 data

library(Seurat)
library(ddqcR)
library(SeuratDisk)

setwd('~/datasets/SDY997')

exprMat <- read.delim('celseq_matrix_ru1_reads.tsv.725701.gz', header=T, row.names = 1)
exprMat[is.na(exprMat)] <- 0
meta <- read.delim('celseq_meta_unfiltered.tsv.725703.gz', header=T)
pbmc <- CreateSeuratObject(counts = exprMat)
pbmc <- subset(pbmc, cells = meta$cell_name)
meta <- subset(meta, cell_name %in% colnames(pbmc))
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

saveRDS(pbmc, 'pbmc.RDS')

# Subset cells for Leukocyte
pbmc <- subset(pbmc, type == 'Leukocyte')

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

cell.count <- table(pbmc$predicted.celltype.l2)
keep <- names(cell.count[cell.count > 100])
pbmc <- subset(pbmc, predicted.celltype.l2 %in% keep)

Idents(pbmc) <- 'predicted.celltype.l2'

pdf('DimPlot.pdf')
DimPlot(pbmc, label = TRUE, reduction='ref.umap', repel=T) + NoLegend()
dev.off()

saveRDS(pbmc, 'pbmc.lymphocytes.RDS')



