library(ggplot2)
library(reshape2)
library(dplyr)

load('../datasets/XCI/chrX.Rdata')
source('../PhD/functions/edgeR.list.R')

pSS <- deg.list('pSS_GSE157278/differential.expression/edgeR/', logfc=0.5)
MS <- deg.list('MS_GSE193770/differential.expression/edgeR/', logfc=0.5)
UC <- deg.list('UC_GSE125527/differential.expression/edgeR/', logfc=0.5)
CD_colon <- deg.list('CD_Kong/colon/differential.expression/edgeR/', logfc=0.5)
CD_TI <- deg.list('CD_Kong/TI/differential.expression/edgeR/', logfc=0.5)
SLE <- deg.list('lupus_Chun/differential.expression/edgeR/', logfc=0.5)

studies <- list(pSS, MS, UC, CD_colon, CD_TI, SLE)

# Heatmap function
heatmap.matrix <- function(study, features=NULL, direction=NULL){
    if(direction == 'up'){
        study <- lapply(study, function(x) subset(x, logFC > 0.5))
    } else if(direction == 'down'){
        study <- lapply(study, function(x) subset(x, logFC < -0.5))
    }
    genes <- unique(Reduce(rbind, study)$gene)
    if(!is.null(features)){
        genes <- genes[genes %in% features]
    }
    celltypes <- names(study)
    mtx <- matrix(NA, nrow=length(genes), ncol=length(celltypes))
    rownames(mtx) <- genes
    colnames(mtx) <- celltypes
    for(cell in celltypes){
        index <- match(study[[cell]]$gene, genes)
        index <- index[!is.na(index)] # remove NAs from index
        if (length(index) > 0){
            mtx[index, cell] <- subset(study[[cell]], gene %in% genes)$logFC
        }
    }
    return(mtx)
}
heatmap.matrix(pSS, rownames(chrX), 'up')
# Build heatmaps for each study
upregulated.matrices <- lapply(studies, function(x) heatmap.matrix(x, rownames(chrX), 'up'))
downregulated.matrices <- lapply(studies, function(x) heatmap.matrix(x, rownames(chrX), 'down'))

# Build heatmaps for each study
pdf('pSS.heatmap.down.pdf')
melt(upregulated.matrices[[1]]) %>%
    ggplot(aes(x=Var1, y=Var2, fill=value) +
    geom_tile() +
    scale_fill_gradient2(low='blue', mid='white', high='red', midpoint=0) +
    labs(title='pSS upregulated') + xlab('') + ylab('')
dev.off()

theme(axis.text.x = element_text('none', family = 'sans'), 
    axis.text.y = element_text(family = 'sans')) +


library(gplots)

pdf('pSS.heatmap.up.pdf')
# Cluster the data by the y-axis
mat <- matrices[[1]]
mat <- mat[rev(order(rowMeans(mat))),]
mat <- mat[,order(colMeans(mat))]

heatmap.2(mat, trace='none', col=rev(redgreen(75)), dendrogram='row',
          key=TRUE, keysize=1.5, key.title='logFC', margins=c(10,10),
          labRow=rownames(mat), labCol=colnames(mat),
          xlab='Cell type', ylab='Gene', main='pSS')
dev.off()