# Load packages
library(Seurat)
library(speckle)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(cluster)
library(ggrepel)
library(fgsea)
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/fishers.test.degs.R')

# Load colours 
load('/directflow/SCCGGroupShare/projects/lacgra/PhD/R/celltype.colours.RData')
load('/directflow/SCCGGroupShare/projects/lacgra/PhD/R/study_colours.RData')

# Create directory to store results
if(!dir.exists('Aim_1_2024')){dir.create('Aim_1_2024')}
if(!dir.exists('Aim_1_2024/figure.data')){dir.create('Aim_1_2024/figure.data')}

# Set variables
pbmc <- readRDS(commandArgs(trailingOnly=TRUE)[1])
disease <- commandArgs(trailingOnly=TRUE)[2]

# Read in gene sets 
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')
disgene <- read.delim(paste0('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/', disease, '.tsv'))$Gene
GWAS <- read.delim(paste0('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/', disease, '_GWAS.tsv'))
GWAS <- unique(unlist(lapply(GWAS$Gene, function(x) unlist(strsplit(x, ';')))))
chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/sex_hormones.RData')

### Figure 2A - UMAP coloured by cell type ###
pdf('Aim_1_2024/Figure_2A.pdf', width = 10, height = 10)
DimPlot(pbmc, group.by='cellTypist', label=TRUE, cols=celltype_colours, order=TRUE)
dev.off()

### Figure 2B - UMAP coloured by condition ###
pdf('Aim_1_2024/Figure_2B.pdf', width = 10, height = 10)
DimPlot(pbmc, group.by='condition', label=FALSE, cols=c('control'='#0000FF', 'disease'='#FF0000'), order=TRUE)
dev.off()

### Figure 3 - Marker gene expression dotplot ###
load('/directflow/SCCGGroupShare/projects/lacgra/PhD/R/marker_genes.RData')
average_expression <- AverageExpression(pbmc, features = marker_genes)$decontXcounts
marker_genes <- rownames(average_expression[rowSums(average_expression) > 0,])

pdf('Aim_1_2024/Figure_3.pdf', width = 20, height = 10)
DotPlot(pbmc, features = marker_genes, cols = c("lightgrey", "red")) + 
theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8)) +
ylab('') + xlab('Marker genes')
dev.off()

### Figure 4A - Cell type proportions clustered bar chart ###
props <- getTransformedProps(clusters = pbmc$cellTypist, 
                             sample = pbmc$individual)

pdf('Aim_1_2024/Figure_4A.pdf', width = 10, height = 10)
p <- plotCellTypeProps(clusters = pbmc$cellTypist, sample = pbmc$individual) + 
    ggtitle("Cell type proportions") + 
    theme(plot.title = element_text(size = 18, hjust = 0)) +
    scale_fill_manual(values = celltype_colours)
p
dev.off()

write.table(p$data, 'Aim_1_2024/figure.data/Figure_4A.txt', sep='\t')

### Figure 4B - Cell type proportions scatter plot ###
pbmc$condition <- factor(pbmc$condition, levels=c('disease', 'control'))
output.logit <- propeller(clusters=pbmc$cellTypist, sample=pbmc$individual, 
    group=pbmc$condition, transform='logit')

# Create a new data frame based on the condition
label_data <- if(nrow(output.logit[output.logit$FDR < 0.05,]) > 0){
    output.logit[output.logit$FDR < 0.05,]
} else {
    output.logit[output.logit$PropRatio > 2,]
}

pdf('Aim_1_2024/Figure_4B.pdf')
ggplot(output.logit, aes(x=PropRatio, y=-log10(FDR), colour=BaselineProp.clusters)) + 
    geom_point() +
    scale_colour_manual(values=celltype_colours) +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    geom_text_repel(data=label_data, aes(label=BaselineProp.clusters), 
                    box.padding = 0.5, point.padding = 0.5, colour='black') +
    theme_minimal() + 
    theme(legend.position = 'none') + 
    xlab('Proportion Ratio') + 
    ylab('-log10(FDR)') + 
    ggtitle('Cell type proportions') + 
    theme(plot.title = element_text(size = 18, hjust = 0))
dev.off()

write.table(output.logit, 'Aim_1_2024/figure.data/Figure_4B.txt', sep='\t')
write.table(label_data, 'Aim_1_2024/figure.data/Figure_4B_label.txt', sep='\t')

### Figure 5A - Barplot of up/downregulated genes ###
degs <- deg.list('differential.expression/edgeR', logfc=0.1)
degs.df <- lapply(degs, function(x){
    up <- sum(x$logFC > 0)
    down <- sum(x$logFC < 0)
    data.frame(up=up, down=down)
})
degs.df <- do.call(rbind, degs.df)
degs.df$down <- degs.df$down * -1
degs.df$celltype <- names(degs)

degs.df <- melt(degs.df, id.vars='celltype')

# Convert variable to factor and reorder levels
degs.df$variable <- factor(degs.df$variable, levels = c("down", "up"))
# Tidy up cell type names 
degs.df$celltype <- replace.names(gsub('_', '.', degs.df$celltype))

pdf('Aim_1_2024/Figure_5A.pdf', width = 10, height = 10)
ggplot(degs.df, aes(x=celltype, y=value, fill=variable)) +
    geom_bar(stat="identity", position="identity") +
    coord_flip() +
    scale_fill_manual(values=c('down'='blue', 'up'='red')) +
    scale_y_continuous(breaks = seq(-100, 200, by = 25)) +
    theme_minimal() +
    labs(x="Cell Type", y="Number of Genes", fill="Direction")
dev.off()

write.table(degs.df, 'Aim_1_2024/figure.data/Figure_5A.txt', sep='\t')

### Figure 5B - Representitive volcano plot ###
top_celltype <- names(sort(unlist(lapply(degs, nrow)), decreasing=TRUE)[1])
example <- read.delim(paste0('differential.expression/edgeR/', top_celltype, '.txt'), sep='\t')

# Select top 10 genes
top_genes <- example %>% 
  arrange(desc(-log10(FDR))) %>% 
    head(10)

pdf('Aim_1_2024/Figure_5B.pdf')
ggplot(example, aes(x=logFC, y=-log10(FDR))) + 
    geom_point(aes(colour=ifelse(FDR < 0.05, 'red', 'black'))) + 
    geom_hline(yintercept=-log10(0.05), linetype='dashed') + 
    geom_vline(xintercept=c(-0.1, 0.1), linetype='dashed') + 
    geom_text_repel(data = top_genes, aes(label = gene), size = 3) +
    theme_minimal() + 
    theme(legend.position = 'none') + 
    xlab('logFC') + 
    ylab('-log10(FDR)') + 
    ggtitle(replace.names(gsub('_', '.', top_celltype))) + 
    theme(plot.title = element_text(size = 18, hjust = 0))
dev.off()

### Figure 6A - Barplot of up/downregulated chrX genes ###
chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name

degs.chrX <- lapply(degs, function(x){
    up.chrX.genes <- subset(x, logFC > 0 & gene %in% chrX)$gene
    up.chrX <- sum(!up.chrX.genes %in% rownames(escape))
    up.escape <- sum(up.chrX.genes %in% rownames(escape))

    down.chrX.genes <- subset(x, logFC < 0 & gene %in% chrX)$gene
    down.chrX <- sum(!down.chrX.genes %in% rownames(escape))
    down.escape <- sum(down.chrX.genes %in% rownames(escape))

    data.frame(up=up.chrX, up.escape=up.escape, down=down.chrX, down.escape=down.escape)
})
degs.chrX <- do.call(rbind, degs.chrX)
degs.chrX$down <- degs.chrX$down * -1
degs.chrX$down.escape <- degs.chrX$down.escape * -1
degs.chrX$celltype <- replace.names(gsub('_', '.', names(degs)))

degs.chrX <- melt(degs.chrX, id.vars='celltype')

# Convert variable to factor and reorder levels
degs.chrX$variable <- factor(degs.chrX$variable, levels = c("down", "down.escape", "up", "up.escape"))

pdf('Aim_1_2024/Figure_6A.pdf')
ggplot(degs.chrX, aes(x=celltype, y=value, fill=variable)) +
    geom_bar(stat="identity", position="stack") +
    coord_flip() +
    scale_y_continuous(labels = abs) +
    theme_minimal() +
    labs(x="Cell Type", y="Number of Genes", fill="Direction") +
    scale_fill_brewer(palette="Set1") +
    ggtitle('Differentially expressed chrX genes')
dev.off()

### Figure 6B - Heatmap of differentially expressed genes across cell types coloured by logFC ###
degs_mtx <- matrix(0, nrow=length(unique(unlist(lapply(degs, function(x) x$gene)))), ncol=length(degs))
rownames(degs_mtx) <- unique(unlist(lapply(degs, function(x) x$gene)))
colnames(degs_mtx) <- replace.names(gsub('_', '.', names(degs)))

for(i in 1:length(degs)){
    degs_mtx[degs[[i]]$gene, i] <- degs[[i]]$logFC
}

pdf('Aim_1_2024/Figure_6B.pdf', width = 8, height = 15)
Heatmap(degs_mtx, name='logFC', col=colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red')), 
        cluster_rows=TRUE, cluster_columns=TRUE, show_row_names=TRUE, row_names_gp = gpar(fontsize = 5) , 
        show_column_names=TRUE, column_names_gp = gpar(fontsize = 12))
dev.off()

### Figure 6C - Fisher's exact test of gene sets ###
edgeR <- deg.list('differential.expression/edgeR', filter=FALSE)

# XCI escape
xcape.enrichment <- lapply(edgeR, function(x){
    tmp <- fisher.test.edgeR(x, logfc = 0.1, rownames(escape))
    tmp$p.value
})
names(xcape.enrichment) <- replace.names(gsub('_', '.', names(xcape.enrichment)))
xcape.enrichment <- data.frame(
    celltype=names(xcape.enrichment), 
    p.value=unlist(xcape.enrichment), 
    FDR=p.adjust(unlist(xcape.enrichment), method='fdr'))

# DisGeNet
disgene.enrichment <- lapply(edgeR, function(x){
    tmp <- fisher.test.edgeR(x, logfc = 0.1, disgene)
    tmp$p.value
})
names(disgene.enrichment) <- replace.names(gsub('_', '.', names(disgene.enrichment)))
disgene.enrichment <- data.frame(
    celltype=names(disgene.enrichment), 
    p.value=unlist(disgene.enrichment),
    FDR=p.adjust(unlist(disgene.enrichment), method='fdr'))

# GWAS
GWAS.enrichment <- lapply(edgeR, function(x){
    tmp <- fisher.test.edgeR(x, logfc = 0.1, GWAS)
    tmp$p.value
})
names(GWAS.enrichment) <- replace.names(gsub('_', '.', names(GWAS.enrichment)))
GWAS.enrichment <- data.frame(
    celltype=names(GWAS.enrichment), 
    p.value=unlist(GWAS.enrichment),
    FDR=p.adjust(unlist(GWAS.enrichment), method='fdr'))

# Estrogen
ER.enrichment <- lapply(edgeR, function(x){
    tmp <- fisher.test.edgeR(x, logfc = 0.1, sex_hormones$ER)
    tmp$p.value
})
names(ER.enrichment) <- replace.names(gsub('_', '.', names(ER.enrichment)))
ER.enrichment <- data.frame(
    celltype=names(ER.enrichment), 
    p.value=unlist(ER.enrichment),
    FDR=p.adjust(unlist(ER.enrichment), method='BH'))

# merge datasets and plot heatmap
combined.enrichment <- data.frame(
    row.names=rownames(xcape.enrichment),
    Escape=xcape.enrichment$FDR, 
    DisGeNet=disgene.enrichment$FDR, 
    GWAS=GWAS.enrichment$FDR, 
    Estrogen=ER.enrichment$FDR)

pdf('Aim_1_2024/Figure_6C.pdf')
Heatmap(as.matrix(combined.enrichment), name='FDR', col=colorRamp2(c(0, 0.06, 1), c('red', 'blue', 'white')), 
        cluster_rows=TRUE, cluster_columns=TRUE, show_row_names=TRUE, row_names_gp = gpar(fontsize = 12), 
        show_column_names=TRUE, column_names_gp = gpar(fontsize = 12))
dev.off()

### GSEA hallmark ###
hallmark <- gmtPathways('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/h.all.v7.5.1.symbols.gmt')

fgsea_list <- list()
for(i in 1:length(edgeR)){
    ranked_genes <- edgeR[[1]]$logFC
    names(ranked_genes) <- edgeR[[i]]$gene

    fgseaRes <- fgsea(pathways = hallmark, 
                    stats    = ranked_genes,
                    minSize  = 15,
                    maxSize  = 500)

    data.frame(subset(fgseaRes, padj < 0.05))$pathway

    collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.05], 
                                        hallmark, ranked_genes)
    mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
                            order(-NES),]
    mainPathways$escape_genes <- unlist(lapply(mainPathways$leadingEdge, function(x) sum(x %in% rownames(escape))))

    mainPathways <- data.frame(mainPathways)

    fgsea_list[[i]] <- mainPathways
}
names(fgsea_list) <- replace.names(gsub('_', '.', names(edgeR)))

save(fgsea_list, file='Aim_1_2024/figure.data/fgsea_list.RData')

fgsea_df <- dplyr::bind_rows(fgsea_list, .id = "celltype")
fgsea_df$celltype <- factor(fgsea_df$celltype)

pdf('Aim_1_2024/Figure_7.pdf')
ggplot(fgsea_df, aes(x=pathway, y=-log10(padj), colour=celltype, size=escape_genes)) + 
    geom_point() + 
    coord_flip() + 
    theme_minimal() +
    scale_color_manual(values=celltype_colours) +
    xlab('Pathway') + 
    ylab('-log10(padj)') + 
    ggtitle('GSEA Hallmark pathways') + 
    theme(plot.title = element_text(size = 18, hjust = 0))
dev.off()

### Figure 9 - Heatmap of SCENIC TF AUC comparison FDR ###
auc <- read.csv('SCENIC/permutation_test_avg.csv')
auc$celltype <- gsub('_.+', '', auc$celltype)
auc$TF <- gsub('.+_', '', auc$TF)

auc_mtx <- reshape2::dcast(auc, TF ~ celltype, value.var='FDR', drop=FALSE, fill=NA)
rownames(auc_mtx) <- auc_mtx$TF

pdf('Aim_1_2024/Figure_9.pdf')
Heatmap(as.matrix(auc_mtx[,-1]), name='FDR', col=colorRamp2(c(0, 0.05, 1), c('red', 'blue', 'white')), 
        cluster_rows=TRUE, cluster_columns=TRUE, show_row_names=TRUE, row_names_gp = gpar(fontsize = 5), 
        show_column_names=TRUE, column_names_gp = gpar(fontsize = 12))
dev.off()