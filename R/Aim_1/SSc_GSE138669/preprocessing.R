library(Seurat)
library(magrittr)
library(ddqcR)
library(reticulate)
library(transformGamPoi)
library(iasva)
library(sva)
library(irlba)
library(SummarizedExperiment)

h5_files <- list.files('h5.files', pattern = '*.h5', full.names = T)

seurat_objects <- lapply(h5_files, Read10X_h5)
names(seurat_objects) <- gsub('GSM\\d+_|raw_feature_bc_matrix.h5', '', basename(h5_files))

# merge the 22 seurat objects
pbmc <- merge(seurat_objects[1], y=seurat_objects[-1], add.cell.ids = names(seurat_objects)[1:2])