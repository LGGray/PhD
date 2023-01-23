# This Rscript reads in a given Seurat object and exports the SCTransformed expression matrix for machine Learning analysis

library(Seurat)

# Take in command line arguments i.e the Seurat object filename
file <- commandArgs(trailingOnly=TRUE)
ifelse(dir.exists('exp.matrix'), print('directory exists'), dir.create('exp.matrix'))

# Read in the Seurat object
pbmc <- readRDS(file)

# Loop through each cell and export the expression matrix
for(cell in levels(pbmc)){
    pbmc.subset <- subset(pbmc, predicted.celltype.l2 == cell)
    class <- pbmc.subset$condition
    exp.matrix <- GetAssayData(pbmc.subset, assay='SCT', slot='counts')
    exp.matrix <- data.frame(t(as.matrix(exp.matrix)))
    exp.matrix <- cbind(class=class, exp.matrix)
    write.table(exp.matrix, paste0('exp.matrix/', gsub(' ', '.', cell), '.txt'), row.names=T, col.names=T, sep='\t', quote=F)
}

print('Complete')
