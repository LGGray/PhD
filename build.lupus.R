library(Seurat)
library(dplyr)
library(ddqcR)

args = commandArgs(trailingOnly=TRUE)
args = as.numeric(args)

setwd('/directflow/SCCGGroupShare/projects/lacgra/lupus.Chun/')

# Rename metadata field
metadata <- read.csv('metadata.csv')
metadata$cell_type %<>% 
  gsub(', ', '_', .) %>%
  gsub(' ', '_', .)
  
# Find input files
files <- list.files(path='expMat', pattern = '.csv', full.names = T)

# Read file and build object
infile <- read.csv(files[args])
rownames(infile) <- infile[,1]
cell <- gsub('.csv', '', basename(files[args]))
pbmc <- CreateSeuratObject(counts=infile[,-1])
pbmc@meta.data <- cbind(pbmc@meta.data, subset(metadata, cell_type %in% cell))


if(all.equal(rownames(pbmc@meta.data), gsub('-', '.', pbmc@meta.data$index))){
  print('Cell names matched')
}

# Remove obvious bad quality cells
pbmc <- initialQC(pbmc)

# Return dataframe of filtering statistics
df.qc <- ddqc.metrics(pbmc)

# Filter out the cells
pbmc <- filterData(pbmc, df.qc)

# SCTransform
pbmc <- SCTransform(pbmc, verbose = FALSE)

#Standard seurat pipeline
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)
pbmc <- FindNeighbors(object = pbmc)
pbmc <- RunUMAP(object=pbmc, dims=1:10)

Idents(pbmc) <- 'cell_type'
saveRDS(pbmc, paste0('seurat.objects/all/', cell, '.RDS'))

female.pbmc <- subset(pbmc, sex == 'female')
saveRDS(female.pbmc, paste0('seurat.objects/female/', cell, '.RDS'))

        
        


