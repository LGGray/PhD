library(Seurat)

# Convert Seurat object into a loom file for SCENIC

pbmc <- readRDS(commandArgs(trailingOnly=TRUE))

pbmc.loom <- as.loom(pbmc, filename = "SCENIC/pbmc.female.loom", verbose = TRUE)
pbmc.loom
pbmc.loom$close_all()