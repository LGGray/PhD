suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scMerge)
  library(scater)
  library(Seurat)
  library(SeuratDisk)
})

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/scMerge')

load('~/datasets/XCI/chrX.Rdata')

# Preprocessing of scRNA data
# 1. ddqcR to remove low quality cells
# 2. Annotate with Azimuth and scPred - find overlap between the two
# 3. Ensure there is are metadata columns named 'individual' and 'condition' to allow overlap between studies
# 4. Set the RNA assay as default assay
# 5. Save objects including both sex and just females i.e pbmc.female.RDS

# Convert seurat objects to SingleCellExperiment
pSS <- as.SingleCellExperiment(readRDS('../pSS_GSE157278/pbmc.female.RDS'))
MS <- as.SingleCellExperiment(readRDS('../MS_GSE193770/pbmc.female.RDS'))
MS.2 <- as.SingleCellExperiment(readRDS('../MS_GSE180759/pbmc.female.RDS'))
RA <- as.SingleCellExperiment(readRDS("../RA_SDY998/pbmc.female.RDS"))
SLE <- as.SingleCellExperiment(readRDS("../SLE_SDY997/pbmc.female.RDS"))
UC <- as.SingleCellExperiment(readRDS("../UC_GSE125527/pbmc.female.RDS"))
UC.2 <- as.SingleCellExperiment(readRDS("../UC_GSE182270/pbmc.female.RDS"))
CD <- as.SingleCellExperiment(readRDS("../CD_GSE134809/pbmc.female.RDS"))
AD <- as.SingleCellExperiment(readRDS("../AD_GSE147424/pbmc.female.RDS"))
# T1D <- as.SingleCellExperiment(readRDS("T1D_parsebio/pbmc.female.RDS"))

# Create list of objects
sce_list = list(
  pSS = pSS,
  MS = MS,
  MS.2 = MS.2,
  RA = RA,
  SLE = SLE,
  UC = UC,
  UC.2 = UC.2,
  CD = CD,
  AD = AD
)

# select column names to keep
keep <- Reduce(intersect, list(names(colData(sce_list[[1]])), names(colData(sce_list[[2]])), names(colData(sce_list[[3]])),
                       names(colData(sce_list[[4]])), names(colData(sce_list[[5]])), names(colData(sce_list[[6]])),
                       names(colData(sce_list[[7]])), names(colData(sce_list[[8]])), names(colData(sce_list[[9]]))))

sce_combine = scMerge::sce_cbind(sce_list = sce_list,
                                 method = "intersect",
                                 colData_names = keep,
                                 batch_names = c("pSS", "MS", "MS.2", "RA", "SLE",
                                                 "UC", "UC.2", "CD", "AD"))
print(sce_combine)

sce_combine = runUMAP(sce_combine, exprs_values = "logcounts")
pdf('logcounts.UMAP.pdf')
scater::plotUMAP(
  sce_combine,
  colour_by = "batch")
dev.off()

# Identify stably expressed genes
exprs_mat = SummarizedExperiment::assay(sce_combine, 'logcounts')
result = scSEGIndex(exprs_mat = exprs_mat, cell_type = sce_combine$predicted.celltype.l2)
# load scMerge SEGs
data(segList)
ctl <- unique(c(segList$human$human_scSEG, rownames(chrX)))


# Identify and remove duplicate cell Ids
dups <- colnames(sce_combine)[which(duplicated(colnames(sce_combine)))]
sce_combine <- sce_combine[,(!colnames(sce_combine) %in% dups)]

scMerge_supervised <- scMerge(
  sce_combine = sce_combine,
  ctl = ctl,
  cell_type = sce_combine$predicted.celltype.l2,
  hvg_exprs = 'logcounts',
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

# SCTransform
merged.seurat <- SCTransform(merged.seurat, verbose = FALSE, assay='originalexp')

Idents(merged.seurat) <- 'predicted.celltype.l2'

pdf('scMerge.DimPlot.pdf')
DimPlot(merged.seurat, label = TRUE, reduction='UMAP', repel=T)
dev.off()

saveRDS(merged.seurat, 'scMerge.RDS')