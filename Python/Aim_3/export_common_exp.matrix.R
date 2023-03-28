library(Seurat)

file <- commandArgs(trailingOnly=TRUE)

pbmc <- readRDS(file)

# Export the expression matrix subsetted by differentially expressed genes for each cell type. Check that there are at least 10 cells per condition
cells <- c('Regulatory T cells', 'Tem/Trm cytotoxic T cells', 'Tcm/Naive helper T cells', 'Tem/effector helper T cells')
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
    exp.matrix <- GetAssayData(pbmc.subset, assay='SCT', slot='counts')
    exp.matrix <- data.frame(t(as.matrix(exp.matrix)))
    # Calculate the correlation matrix to identify dependent features
    cor_matrix <- cor(exp.matrix, method='spearman')
    # identify pairs of features with correlation coefficient > 0.7
    high_cor <- which(abs(cor_matrix) > 0.7 & upper.tri(cor_matrix), arr.ind = TRUE)
    # remove the features with high correlation coefficients
    exp.matrix <- if(nrow(high_cor) > 0){
        print(paste('Removing', nrow(high_cor), 'highly correlated features from', cell, sep=' '))
        exp.matrix <- exp.matrix [,-unique(high_cor[,2]),]
        exp.matrix
    } else {
        print(paste('No highly correlated features in', cell, sep=' '))
        exp.matrix
    }
    exp.matrix <- cbind(class=class, exp.matrix)
    cell <- gsub('/| |-', '.', cell)
    saveRDS(exp.matrix, paste0('exp.matrix/', cell, '.common.RDS'))
}