library(ggplot2)
library(dplyr)
library(Seurat)
library(car)
library(rstatix)

if(dir.exists('variance') != TRUE){dir.create('variance')}

pbmc <- readRDS('pbmc.female.RDS')

result <- lapply(levels(pbmc), function(cell){
    # Select cell type
    print(cell)
    # subset object by cell type
    pbmc.cell <- subset(pbmc, cellTypist == cell)
    # check if there are enough cell in both conditions and skip if not
    if(length(unique(pbmc.cell$condition)) != 2){
        return("Not enough conditions")
        next
    }
    # exp <- GetAssayData(pbmc.cell)
    # keep <- apply(exp, 1, function(x) sum(x > 10) > ncol(pbmc)*0.05) 
    # features <- names(keep[keep == T])
    # pbmc.cell <- subset(pbmc.cell, features=features)

    # extract expression for both conditions
    control <- GetAssayData(subset(pbmc.cell, condition=='control'), slot='counts')
    disease <- GetAssayData(subset(pbmc.cell, condition=='disease'), slot='counts')
    # identify genes with expression in both condition
    control.features <- names(which(rowSums(control) > 0))
    disease.features <- names(which(rowSums(disease) > 0))
    features <- intersect(control.features, disease.features)
    # Extract new features
    control <- GetAssayData(subset(pbmc.cell, condition=='control', features=features), slot='counts')
    disease <- GetAssayData(subset(pbmc.cell, condition=='disease', features=features), slot='counts')

    variance.test <- lapply(1:nrow(control), function(g){
        df <- data.frame(condition=c(rep('control', ncol(control)), rep('disease', ncol(disease))), 
                     expr=c(control[g,], disease[g,]))
        oneway.test(expr ~ condition, data=df, var.equal = FALSE)
        })
    names(variance.test) <- rownames(control)

    tmp <- lapply(names(variance.test), function(gene_name){
        x <- variance.test[[gene_name]]
        data.frame(
        gene=gene_name, 
        F=x$statistic, 
        pvalue=x$p.value)
    })
    tmp <- dplyr::bind_rows(tmp)
    tmp$FDR <- p.adjust(tmp$pvalue, method='fdr')
    write.table(tmp, file=paste0('variance/', gsub(' |/|-', '_', cell), '.txt'), sep='\t', row.names=F, quote=F)
    return(tmp)
})
names(result) <- levels(pbmc)

