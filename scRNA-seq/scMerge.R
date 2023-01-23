suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scMerge)
  library(scater)
  library(Seurat)
})

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/scMerge')

# Preprocessing of scRNA data
# 1. ddqcR to remove low quality cells
# 2. Annotate with Azimuth and scPred - find overlap between the two
# 3. Ensure there is a metadata column named 'condition' to allow overlap between studies
# 4. Set the RNA assay as default assay
# 5. Save objects including both sex and just females i.e pbmc.female.RDS

# Convert seurat objects to SingleCellExperiment
pSS <- as.SingleCellExperiment(readRDS('pSS_GSE157278/pbmc.RDS'))
MS <- as.SingleCellExperiment(readRDS('MS_GSE193770/pbmc.female.RDS'))
RA <- as.SingleCellExperiment(readRDS("RA_SDY998/pbmc.female.RDS"))
SLE <- as.SingleCellExperiment(readRDS("SLE_SDY997/pbmc.lymphocytes.female.RDS"))
UC <- as.SingleCellExperiment(readRDS("UC_GSE125527/pbmc.female.RDS"))
CD <- as.SingleCellExperiment(readRDS("CD_GSE134809/pbmc.female.RDS"))
T1D <- as.SingleCellExperiment(readRDS("T1D_parsebio/pbmc.female.RDS"))

# Create list of objects
sce_list = list(
  pSS = pSS,
  MS = MS,
  RA = RA,
  SLE = SLE,
  UC = UC,
  CD = CD,
  T1D = T1D
)

# select column names to keep
keep <- intersect(names(colData(sce_list[[1]])), names(colData(sce_list[[2]])))

sce_combine = scMerge::sce_cbind(sce_list = sce_list,
                                 method = "intersect",
                                 colData_names = keep,
                                 batch_names = c("pSS", "MS", "RA", "SLE",
                                                 "UC", "CD", "T1D"))
sce_combine

sce_combine = runUMAP(sce_combine, exprs_values = "logcounts")
pdf('logcounts.UMAP.pdf')
scater::plotUMAP(
  sce_combine,
  colour_by = "batch")
dev.off()

# Identify stably expressed genes
exprs_mat = SummarizedExperiment::assay(sce_combine, 'counts')
result = scSEGIndex(exprs_mat = exprs_mat)

scMerge_supervised <- scMerge(
  sce_combine = sce_combine,
  ctl = rownames(result),
  cell_type = sce_combine$predicted.celltype.l2,
  replicate_prop = 1,
  assay_name = "scMerge_supervised",
  verbose = TRUE
)

scMerge_supervised = runUMAP(scMerge_supervised, exprs_values = "scMerge_supervised")

pdf('scMerge.batch.UMAP.pdf')
scater::plotUMAP(
  scMerge_supervised,
  colour_by = "batch")
dev.off()

pdf('scMerge.UMAP.celltype.pdf')
scater::plotUMAP(
  scMerge_supervised,
  colour_by = "predicted.celltype.l2")
dev.off()

# Convert SCE to Seurat
merged.seurat <- as.Seurat(scMerge_supervised)
saveRDS(merged.seurat, 'scMerge.RDS')
