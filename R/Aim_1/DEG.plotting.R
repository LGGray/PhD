library(dplyr)
library(ggplot2)
library(reshape2)
library(gplots)

source('../../../PhD/functions/edgeR.list.R')
source('../../../PhD/functions/chisq.test.degs.R')

load('../../../datasets/XCI/chrX.Rdata')
load('../../../datasets/XCI/escapees.Rdata')
Xi_drivers <- read.delim('../../../datasets/XCI/Xi_drivers.txt', header=F, sep='\t')

edgeR <- deg.list('edgeR/', filter=F)

chrX.enrichment <- lapply(edgeR, function(x){
    tmp <- chisq.test.edgeR(x, rownames(chrX), 0.05)
    data.frame(p=tmp$p.value, size=nrow(subset(x, abs(logFC) > 0.05 & FDR < 0.05 & gene %in% rownames(chrX))))
})
chrX.enrichment <- dplyr::bind_rows(chrX.enrichment, .id='celltype')
chrX.enrichment$FDR <- p.adjust(chrX.enrichment$p, method='fdr')
subset(chrX.enrichment, FDR < 0.05)

pdf('chrX_enrichment.pdf', width=10, height=10)
ggplot(chrX.enrichment, aes(x=celltype, y=-log10(p))) + geom_point(aes(size=size, color=-log10(p))) +
    scale_color_gradient(low="white", high="red") + theme_bw() + theme(panel.grid.major = element_blank()) +
    scale_size(name = "# chrX") +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + ylab("-log10(p)") + xlab("Cell type") + ggtitle("chrX enrichment") +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

escape.enrichment <- lapply(edgeR, function(x){
    tmp <- chisq.test.edgeR(x, rownames(escape), 0.5)
    data.frame(p=tmp$p.value, size=nrow(subset(x, abs(logFC) > 0.5 & FDR < 0.05 & gene %in% rownames(escape))))
})
escape.enrichment <- bind_rows(escape.enrichment, .id='celltype')
pdf('escape_enrichment.pdf', width=10, height=10)
ggplot(escape.enrichment, aes(x=celltype, y=-log10(p))) + geom_point(aes(size=size, color=-log10(p))) +
    scale_color_gradient(low="white", high="red") + theme_bw() + theme(panel.grid.major = element_blank()) +
    scale_size(name = "# escapees") +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + ylab("-log10(p)") + xlab("Cell type") + ggtitle("escapee enrichment") +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Creating heatmap of chrX across cell types
deg.matrix <- matrix(0, nrow=length(rownames(chrX)), ncol=length(edgeR))
rownames(deg.matrix) <- rownames(chrX)
colnames(deg.matrix) <- names(edgeR)
for (i in 1:length(edgeR)){
    tmp <- subset(edgeR[[i]], abs(logFC) > 0.05 & FDR < 0.05 & gene %in% rownames(chrX))
    # match tmp$genes to deg.matrix and fill deg.matrix[,i]] with tmp$logFC
    deg.matrix[match(tmp$gene, rownames(chrX)), i] <- tmp$logFC
}
# filter out genes that are not differentially expressed in any cell type
deg.matrix <- deg.matrix[rowSums(deg.matrix) != 0,]

# plot heatmap
pdf('chrX_heatmap.pdf', width=10, height=10)
heatmap.2(deg.matrix, trace='none', col=rev(colorRampPalette(c("blue", "white", "red"))(100)), 
          scale='none', margins=c(10,10), main='Differentially expressed chrX', xlab='', ylab='',
          srtCol=45, key=FALSE)
dev.off()

edgeR.deg <- deg.list('edgeR/', logfc=0.5)
deg.df <- lapply(edgeR.deg, function(x){
    data.frame(upregulated=sum(x$logFC > 0), downregulated=sum(x$logFC < 0))
})
deg.df <- bind_rows(deg.df, .id='celltype')

# Plot colum plot of up/down regulated genes
pdf('DEG.barplot.pdf', width=10, height=10)
ggplot(melt(deg.df), aes(x=celltype, y=value, fill=variable)) +
    geom_col(position='dodge') +
    scale_fill_manual(values=c('upregulated'='red', 'downregulated'='blue')) +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + ylab("# DEGs") +
    ggtitle("Differentially expressed genes") +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

library(clusterProfiler)

edgeR <- deg.list('pSS_GSE157278/differential.expression/edgeR/', filter=F)
gene.set <- read.gmt('../gene.sets/h.all.v7.5.1.symbols.gmt')
ORA <- enricher(gene = subset(edgeR[[1]], abs(logFC) > 0.5 & FDR < 0.05)$gene, 
                   universe = edgeR[[1]]$gene,
                   pAdjustMethod = "fdr", 
                   qvalueCutoff = 0.05,
                   TERM2GENE = gene.set)
result <- ORA@result
result <- subset(result, qvalue < 0.05)
keep <- lapply(result$geneID, function(x){
    any(rownames(chrX) %in% unlist(strsplit(x, '/')))
})
result.chrX <- result[unlist(keep),]
genes <- unlist(strsplit(result.chrX$geneID, '/'))
genes[genes %in% rownames(chrX)]
