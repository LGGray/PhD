library(Seurat)

file <- commandArgs(trailingOnly=TRUE)

pbmc <- readRDS(file)

# Export the expression matrix subsetted by differentially expressed genes for each cell type. Check that there are at least 10 cells per condition
cells <- c('Regulatory T cells', 'Tem/Trm cytotoxic T cells', 'Tcm/Naive helper T cells', 'Tem/Effector helper T cells')
for(cell in cells){
    if(cell %in% levels(pbmc) == FALSE){
        cat('No cells of this type. Skipping', cell, '\n')
        next
    }
    common <- read.delim(paste0('../../common.features/', gsub(' |-|/', '.', cell), '.chrX.txt'), header=F)
    pbmc.subset <- subset(pbmc, cellTypist == cell, features=common$V1)
    if (length(which(pbmc.subset$condition == 'control')) < 10 | length(which(pbmc.subset$condition == 'disease')) < 10) {
        cat('Not enough samples. Skipping', cell, '\n')
        next
    }
    class <- pbmc.subset$condition
    individual <- pbmc.subset$individual
    exp.matrix <- GetAssayData(pbmc.subset, assay='SCT', slot='counts')
    exp.matrix <- data.frame(t(as.matrix(exp.matrix)))
    # Calculate the correlation matrix to identify dependent features
    cor_matrix <- cor(exp.matrix, method='spearman')
    # identify pairs of features with correlation coefficient > 0.9
    high_cor <- which(abs(cor_matrix) > 0.9 & upper.tri(cor_matrix), arr.ind = TRUE)
    # remove the features with high correlation coefficients
    exp.matrix <- if(nrow(high_cor) > 0){
        features <- rownames(high_cor)
        exp.matrix.sub <- exp.matrix[,features]
        pca <- prcomp(exp.matrix.sub, scale=TRUE)
        exp.matrix <- exp.matrix[,!(colnames(exp.matrix) %in% features)]
        exp.matrix[, paste(features, collapse='_')] <- pca$x[,1]
        print(paste('Reducing dimensions of', nrow(high_cor), 'highly correlated features from', cell, sep=' '))
    } else {
        print(paste('No highly correlated features in', cell, sep=' '))
        exp.matrix
    }
    exp.matrix <- cbind(class=class, individual=individual, exp.matrix)
    cell <- gsub('/| |-', '.', cell)
    saveRDS(exp.matrix, paste0('exp.matrix/', cell, '.common.RDS'))
}