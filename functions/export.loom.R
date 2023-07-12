library(Seurat)
library(SeuratDisk)

pbmc <- readRDS('pbmc.female.RDS')

pbmc.loom <- as.loom(pbmc, filename = "pbmc.loom", verbose = FALSE)
pbmc.loom
pbmc.loom$close_all()

