library(Seurat)
library(SeuratDisk)

setwd('~/datasets/lupus.Chun')

disease <- LoadH5Seurat('disease.h5seurat')
Idents(disease) <- 'cell_type'
control <- LoadH5Seurat('control.h5seurat')
Idents(control) <- 'cell_type'

# 1-11

args = commandArgs(trailingOnly=TRUE)
args <- as.numeric(args)

dis <- subset(disease, cell_type %in% levels(disease)[args])
cnt <- subset(control, cell_type %in% levels(disease)[args])
pbmc <- merge(dis, cnt)
saveRDS(pbmc, paste0('seurat.object/',
                     gsub(' |, ', '_', levels(disease)[args]), '.RDS'))

  