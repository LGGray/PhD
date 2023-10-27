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


##############
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(ggplot2)

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')

edgeR <- deg.list('differential.expression/edgeR', filter=F)
deg <- deg.list('differential.expression/edgeR', logfc=0.5)

deg.escape <- lapply(deg, function(x) x[x$gene %in% rownames(escape),])

# Create matrix of logFC.disease_vs_control for each gene
genes <- unique(unlist(lapply(deg.escape, function(x) x$gene)))
mtx <- matrix(0, nrow=length(genes), ncol=length(deg.escape))
rownames(mtx) <- genes
colnames(mtx) <- replace.names(gsub('_', '.', names(deg.escape)))
for (i in 1:length(deg.escape)) {
  mtx[,i] <- deg.escape[[i]]$logFC.disease_vs_control[match(genes, deg.escape[[i]]$gene)]
}
mtx[is.na(mtx)] <- 0

mtx.cor <- cor(mtx, method='spearman')
mtx.cor[is.na(mtx.cor)] <- 0

# Cor heatmap
pdf('APR/escape.cor.pdf')
Heatmap(mtx.cor, clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
clustering_method_rows = "complete", clustering_method_columns = "complete", column_title = 'Primary Sjogren\'s Syndrome',
col=colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), name="Spearman's rho", column_title_side = "top",
row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
dev.off()

# PCA logFC matrix
pca <- prcomp(t(mtx), center=T)
pdf('APR/chrX.pca.pdf')
ggplot(data.frame(pca$x), aes(x=PC1, y=PC2)) + geom_point() + theme_bw() +
geom_text(aes(label=rownames(pca$x)), size=2, hjust=0, vjust=0)
dev.off()

library(biomaRt)

rownames(mtx) <- gsub('FAM122B', 'PABIR2', rownames(mtx))
rownames(mtx) <- gsub('CXorf40A', 'EOLA1', rownames(mtx))
rownames(mtx) <- gsub('CXorf40B', 'EOLA2', rownames(mtx))
rownames(mtx) <- gsub('TAZ', 'TAFAZZIN', rownames(mtx))

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene_info <- getBM(attributes = c("start_position", "external_gene_name", "chromosome_name"), filters = "external_gene_name", values = rownames(mtx), mart = mart)
gene_info <- subset(gene_info, chromosome_name == 'X')
gene_info <- gene_info[order(gene_info$start_position, decreasing=FALSE),]

mtx.ordered <- mtx[match(gene_info$external_gene_name, rownames(mtx)),]

pdf('APR/chrX.heatmap.pdf')
Heatmap(mtx.ordered, cluster_rows = FALSE, name="logFC", column_title = '',
col=colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")), rect_gp = gpar(col = "white", lwd = 2),
show_row_names = FALSE)
dev.off()
