library(Seurat)
library(SeuratDisk)
library(data.table)

pbmc <- readRDS(commandArgs(trailingOnly=TRUE))

# Export counts as .csv
expr_matrix <- as.matrix(GetAssayData(pbmc, slot = 'data'))
fwrite(data.frame(t(expr_matrix)), 'SCENIC/expr_matrix.csv', row.names = TRUE)

# Export counts as loom
pbmc.loom <- as.loom(pbmc, filename = "SCENIC/pbmc.female.loom", verbose = TRUE, overwrite = TRUE)
pbmc.loom
pbmc.loom$close_all()

SaveH5Seurat(pbmc, filename = "pbmc.female.control-managed.h5Seurat")
Convert("pbmc.female.control-managed.h5Seurat", dest = "h5ad")