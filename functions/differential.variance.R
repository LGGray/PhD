library(ggplot2)
library(dplyr)
library(Seurat)
library(car)
library(rstatix)

source('../../PhD/functions/chisq.test.degs.R')
load('../../datasets/XCI/chrX.Rdata')

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
    exp <- GetAssayData(pbmc.cell)
    keep <- apply(exp, 1, function(x) sum(x > 0) > ncol(pbmc)*0.05) 
    features <- names(keep[keep == T])
    pbmc.cell <- subset(pbmc.cell, features=features)

    control <- GetAssayData(subset(pbmc.cell, condition=='control'), slot='counts')
    disease <- GetAssayData(subset(pbmc.cell, condition=='disease'), slot='counts')

    cl <- makeCluster(detectCores()-50)
    clusterExport(cl, c("control", "disease", "oneway.test", "games_howell_test"))
    variance.test <- parLapply(cl, 1:nrow(control), function(g){
        df <- data.frame(condition=c(rep('control', ncol(control)), rep('disease', ncol(disease))), 
                     expr=c(control[g,], disease[g,]))
        oneway.test(expr ~ condition, data=df, var.equal = FALSE)
        # games_howell_test(expr ~ condition, data=df)
    })
    names(variance.test) <- rownames(control)
    stopCluster(cl)

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

enrichment.test <- function(data, genes)
{
  if (is.data.frame(data) == F) stop("data is not a data frame")
  if (is.vector(genes) == F) stop("genes is not a vector")
  a <- length(intersect(data[data$FDR < 0.05,'gene'], genes))
  b <- length(setdiff(data[data$FDR < 0.05,'gene'], genes))
  c <- length(intersect(data[data$FDR > 0.05,'gene'], genes))
  d <- length(setdiff(data[data$FDR > 0.05,'gene'], genes))
  test = (chisq.test(matrix(c(a,b,c,d), nrow=2)))
  return(test)
}
print('Enrichment test')
lapply(result, function(x) enrichment.test(x, rownames(chrX))$p.value)

