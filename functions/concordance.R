# Code to determine concordance of X-linked gene expression between studies

source('../PhD/functions/edgeR.list.R')
load('../datasets/XCI/chrX.Rdata')

pSS <- deg.list('pSS_GSE157278/differential.expression/edgeR', logfc=0.05)
UC <- deg.list('UC_GSE125527/differential.expression/edgeR', logfc=0.05)
CD_colon <- deg.list('CD_Kong/colon/differential.expression/edgeR', logfc=0.05)
CD_TI <- deg.list('CD_Kong/TI/differential.expression/edgeR', logfc=0.05)
SLE <- deg.list('lupus_Chun/differential.expression/edgeR', logfc=0.05)
studies <- list('pSS'=pSS, 'UC'=UC, 'CD_colon'=CD_colon, 'CD_TI'=CD_TI, 'SLE'=SLE)

# Common celltypes
celltypes <- Reduce(intersect, list(names(pSS), names(UC), names(CD_colon), names(CD_TI), names(SLE)))

# Create matrix of logFC for each celltype
matrix.list <- list()
for(cell in celltypes){
    genes <- unique(unlist(lapply(studies, function(x) subset(x[[cell]], gene %in% rownames(chrX))$gene)))
    plot.matrix <- matrix(0, nrow=length(genes), ncol=length(studies))
    rownames(plot.matrix) <- genes
    colnames(plot.matrix) <- names(studies)
    # Match genes to rownames
    for (i in 1:length(studies)){
        df <- subset(studies[[i]][[cell]], gene %in% rownames(chrX))
        plot.matrix[match(df$gene, genes),i] <- df$logFC.disease_vs_control
    }
    matrix.list[[cell]] <- plot.matrix
}

correlation <- lapply(matrix.list, function(x) cor(x, method='spearman'))

# Heatmap of correlation
library(gplots)
pdf('Tem_Trm_cytotoxic_T_cells.cor.heatmap.pdf')
heatmap.2(correlation[[6]], trace='none', col=colorRampPalette(c('blue', 'white', 'red'))(100), 
          scale='none', dendrogram='none', key=F, cexRow=0.5, srtCol=45, cexCol=0.5,  margins=c(10,10))
dev.off()

# Jaccard distance within studies
library(proxy)
# loop over the cell types in matrix.list
for (i in 1:length(matrix.list)) {
  # get the matrix of logFC values for the i-th cell type
  
  # convert the logFC values to 0 or 1
  binary_matrix <- ifelse(matrix.list[[i]] == 0, 0, 1)
  
  # compute the Jaccard distance between studies
  jaccard_sim <- proxy::simil(binary_matrix, by_rows = FALSE, method = "Jaccard")
  
  # store the Jaccard distance matrix in jaccard_list with the same name as the cell type
  jaccard_list[[names(matrix.list)[i]]] <- jaccard_dist
}


# Correlation within studies
matrix.list <- list()
for(study in studies){
    genes <- unique(unlist(lapply(study, function(x) subset(x, gene %in% rownames(chrX))$gene)))
    plot.matrix <- matrix(0, nrow=length(genes), ncol=length(study))
    rownames(plot.matrix) <- genes
    colnames(plot.matrix) <- names(study)
    # Match genes to rownames
    for (i in 1:length(study)){
        df <- subset(study[[i]], gene %in% rownames(chrX))
        plot.matrix[match(df$gene, genes),i] <- df$logFC.disease_vs_control
    }
    correlation <- cor(plot.matrix, method='spearman')
    correlation[is.na(correlation)] <- 0
    matrix.list <- append(matrix.list, list(correlation))
}

# Heatmap of correlation
pdf('SLE.cor.heatmap.pdf')
heatmap.2(matrix.list[[5]], trace='none', col=colorRampPalette(c('blue', 'white', 'red'))(100), 
          scale='none', dendrogram='none', key=F, cexRow=0.5, srtCol=45, cexCol=0.5,  margins=c(10,10))
dev.off()

# Jaccard distance of genes between studies

