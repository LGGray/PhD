# This Rscript reads in a given Seurat object and exports the SCTransformed expression matrix for machine Learning analysis
options(warn=-1)
library(Seurat)

# Load X chromosome genes
load('../../datasets/XCI/chrX.Rdata')


# Take in command line arguments i.e the Seurat object filename
file <- commandArgs(trailingOnly=TRUE)
directory <- dirname(file)

# Set the working directory
setwd(directory)

# Create the directory to store the expression matrices
ifelse(dir.exists('exp.matrix'), 'directory exists', dir.create('exp.matrix'))

# Read in the Seurat object
pbmc <- readRDS(file)

# Tabulate the number of cells per cell type for each condition
print(table(pbmc$condition, pbmc$cellTypist))

# Export the complete expression matrix for each cell type. Check that there are at least 10 cells per condition
for(cell in levels(pbmc)){
    pbmc.subset <- subset(pbmc, cellTypist == cell)
    if (length(which(pbmc.subset$condition == 'control')) < 10 | length(which(pbmc.subset$condition == 'disease')) < 10) {
        cat('Not enough samples. Skipping', cell, '\n')
        next
    }
    class <- pbmc.subset$condition
    exp.matrix <- GetAssayData(pbmc.subset, assay='SCT', slot='counts')
    exp.matrix <- data.frame(t(as.matrix(exp.matrix)))
    exp.matrix <- cbind(class=class, exp.matrix)
    cell <- gsub('/| |-', '.', cell)
    saveRDS(exp.matrix, paste0('exp.matrix/', cell, '.RDS'))
}

print('Complete Matrix Exported')

# Export the expression matrix subsetted by differentially expressed genes for each cell type. Check that there are at least 10 cells per condition
for(cell in levels(pbmc)){
    files = list.files('psuedobulk', pattern = 'edgeR-LRT.txt')
    if(paste0(gsub('/|-| ', '_', cell), '.edgeR-LRT.txt') %in% files == FALSE){
        cat('No DEGs. Skipping', cell)
        next
    }
    deg <- read.delim(paste0('psuedobulk/', gsub('/|-| ', '_', cell), '.edgeR-LRT.txt'))
    genes <- subset(deg, FDR < 0.05 & abs(logFC) > 0.5)$gene
    pbmc.subset <- subset(pbmc, cellTypist == cell, features=genes)
    if (length(which(pbmc.subset$condition == 'control')) < 10 | length(which(pbmc.subset$condition == 'disease')) < 10) {
        cat('Not enough samples. Skipping', cell, '\n')
        next
    }
    class <- pbmc.subset$condition
    exp.matrix <- GetAssayData(pbmc.subset, assay='SCT', slot='counts')
    exp.matrix <- data.frame(t(as.matrix(exp.matrix)))
    exp.matrix <- cbind(class=class, exp.matrix)
    cell <- gsub('/| |-', '.', cell)
    saveRDS(exp.matrix, paste0('exp.matrix/', cell, '.deg.RDS'))
}

print('DEG Matrix Exported')

# Export the expression matrix subsetted by X chromosome genes for each cell type. Check that there are at least 10 cells per condition
for(cell in levels(pbmc)){
    pbmc.subset <- subset(pbmc, cellTypist == cell, features=rownames(chrX))
    if (length(which(pbmc.subset$condition == 'control')) < 10 | length(which(pbmc.subset$condition == 'disease')) < 10) {
        cat('Not enough samples. Skipping', cell, '\n')
        next
    }
    class <- pbmc.subset$condition
    exp.matrix <- GetAssayData(pbmc.subset, assay='SCT', slot='counts')
    exp.matrix <- data.frame(t(as.matrix(exp.matrix)))
    exp.matrix <- cbind(class=class, exp.matrix)
    cell <- gsub('/| |-', '.', cell)
    saveRDS(exp.matrix, paste0('exp.matrix/', cell, '.chrX.RDS'))
}

print('chrX Matrix Exported')