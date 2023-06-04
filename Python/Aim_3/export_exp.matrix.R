# This Rscript reads in a given Seurat object and exports the SCTransformed expression matrix for machine Learning analysis
options(warn=-1)
library(Seurat)

# Load X chromosome genes
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')

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

# Export the expression matrix subsetted by X chromosome genes for each cell type. Check that there are at least 10 cells per condition
min_cells_per_condition <- 10
for(cell in levels(pbmc)){
    pbmc.subset <- subset(pbmc, cellTypist == cell, features=rownames(chrX))
    if (length(which(pbmc.subset$condition == 'control')) < min_cells_per_condition | length(which(pbmc.subset$condition == 'disease')) < min_cells_per_condition) {
        cat('Not enough samples. Skipping', cell, '\n')
        next
    }
    class <- pbmc.subset$condition
    individual <- pbmc.subset$individual
    exp.matrix <- GetAssayData(pbmc.subset, assay='RNA', slot='counts')
    exp.matrix <- data.frame(t(as.matrix(exp.matrix)))
    # # Calculate the correlation matrix to identify dependent features
    # cor_matrix <- cor(exp.matrix, method='spearman')
    # # identify pairs of features with correlation coefficient > 0.9
    # high_cor <- which(abs(cor_matrix) > 0.9 & upper.tri(cor_matrix), arr.ind = TRUE)
    # # remove the features with high correlation coefficients
    # exp.matrix <- if(nrow(high_cor) > 0){
    #     df <- data.frame(X=rownames(cor_matrix)[high_cor[,1]], Y=colnames(cor_matrix)[high_cor[,2]])
    #     for(x in 1:nrow(df)){
    #         exp.matrix[, paste(df[x,], collapse='_')] <- exp.matrix[, df[x,1]] * exp.matrix[, df[x,2]]
    #     }
    #     features <- unique(unlist(df))
    #     print(paste('Creating interaction term from', nrow(high_cor), 'highly correlated features within', cell, sep=' '))
    #     exp.matrix <- exp.matrix[,!(colnames(exp.matrix) %in% features)]
    #     exp.matrix
    # } else {
    #     print(paste('No highly correlated features in', cell, sep=' '))
    #     exp.matrix
    # }
    exp.matrix <- cbind(class=class, individual=individual, exp.matrix)
    cell <- gsub('/| |-', '.', cell)
    saveRDS(exp.matrix, paste0('exp.matrix/', cell, '.chrX.RDS'))
}

print('chrX Matrix Exported')

# Export the expression matrix subsetted by top 2000 HVG for each cell type. Check that there are at least 10 cells per condition
min_cells_per_condition <- 10
for(cell in levels(pbmc)){
    pbmc.subset <- subset(pbmc, cellTypist == cell)
    if (length(which(pbmc.subset$condition == 'control')) < min_cells_per_condition | length(which(pbmc.subset$condition == 'disease')) < min_cells_per_condition) {
        cat('Not enough samples. Skipping', cell, '\n')
        next
    }
    class <- pbmc.subset$condition
    individual <- pbmc.subset$individual
    pbmc.subset <- FindVariableFeatures(pbmc.subset, nfeatures=2000)
    pbmc.subset <- subset(pbmc.subset, features = VariableFeatures(pbmc.subset))
    exp.matrix <- GetAssayData(pbmc.subset, assay='RNA', slot='counts')
    exp.matrix <- data.frame(t(as.matrix(exp.matrix)))
    exp.matrix <- cbind(class=class, individual=individual, exp.matrix)
    cell <- gsub('/| |-', '.', cell)
    saveRDS(exp.matrix, paste0('exp.matrix/', cell, '.HVG.RDS'))
}

print('HVG Matrix Exported')

# # Create directories to store ML results
ifelse(dir.exists('ML.models'), 'directory exists', dir.create('ML.models'))
ifelse(dir.exists('exp.matrix/AUROC'), 'directory exists', dir.create('exp.matrix/AUROC'))
ifelse(dir.exists('exp.matrix/PRC'), 'directory exists', dir.create('exp.matrix/PRC'))
ifelse(dir.exists('exp.matrix/metrics'), 'directory exists', dir.create('exp.matrix/metrics'))

files <- dir('exp.matrix', pattern='.RDS', full.names=TRUE)
write.table(files, 'exp.matrix/file.list.txt', quote=FALSE, row.names=FALSE, col.names=FALSE)