### Load packages ###
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
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')

### Load colours ###
load('/directflow/SCCGGroupShare/projects/lacgra/PhD/R/celltype.colours.RData')
load('/directflow/SCCGGroupShare/projects/lacgra/PhD/R/study_colours.RData')

### Create directory to store results ###
if(!dir.exists('Aim_1_2024')){dir.create('Aim_1_2024')}

### Set variables ###
pbmc <- readRDS(commandArgs(trailingOnly=TRUE))

# Figure 2A - UMAP coloured by cell type
pdf('Aim_1_2024/Figure_2A.pdf', width = 10, height = 10)
DimPlot(pbmc, group.by='cellTypist', label=TRUE, cols=celltype_colours, order=TRUE)
dev.off()

# Figure 2B - UMAP coloured by condition
pdf('Aim_1_2024/Figure_2B.pdf', width = 10, height = 10)
DimPlot(pbmc, group.by='condition', label=FALSE, cols=c('control'='#0000FF', 'disease'='#FF0000'), order=TRUE)
dev.off()

# Figure 3 - Marker gene expression dotplot
load('/directflow/SCCGGroupShare/projects/lacgra/PhD/R/marker_genes.RData')
average_expression <- AverageExpression(pbmc, features = marker_genes)$decontXcounts
marker_genes <- rownames(average_expression[rowSums(average_expression) > 0,])

pdf('Aim_1_2024/Figure_3.pdf', width = 20, height = 10)
DotPlot(pbmc, features = marker_genes, cols = c("lightgrey", "red")) + 
theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8)) +
ylab('') + xlab('Marker genes')
dev.off()

# Figure 4A - Cell type proportions clustered bar chart
props <- getTransformedProps(clusters = pbmc$cellTypist, 
                             sample = pbmc$individual)

pdf('Aim_1_2024/Figure_4A.pdf', width = 10, height = 10)
plotCellTypeProps(clusters = pbmc$cellTypist, sample = pbmc$individual) + 
    ggtitle("Cell type proportions") + 
    theme(plot.title = element_text(size = 18, hjust = 0)) +
    scale_fill_manual(values = celltype_colours)
dev.off()

# Figure 4B - Cell type proportions scatter plot
pbmc$condition <- factor(pbmc$condition, levels=c('disease', 'control'))
output.logit <- propeller(clusters=pbmc$cellTypist, sample=pbmc$individual, 
    group=pbmc$condition, transform='logit')

# Create a new data frame based on the condition
label_data <- if(nrow(output.logit[output.logit$FDR < 0.05]) > 0){
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

# Figure 5A - Barplot of up/downregulated genes
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
    scale_y_continuous(labels = abs) +
    theme_minimal() +
    labs(x="Cell Type", y="Number of Genes", fill="Direction") +
    scale_fill_brewer(palette="Set1")
dev.off()

# Figure 5B - Representitive volcano plot

top_celltype <- names(sort(unlist(lapply(degs, nrow)), decreasing=TRUE)[1])

example <- read.delim(paste0('differential.expression/edgeR/', top_celltype, '.txt'), sep='\t')

library(ggrepel)

# Select top 10 genes
top_genes <- example %>% 
  arrange(desc(-log10(FDR))) %>% 
    head(10)
    
top_genes <- example[order(example$FDR, decreasing=FALSE),][1:10,'gene']

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

# Figure 6A - Barplot of up/downregulated chrX genes
chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name

degs.chrX <- lapply(degs, function(x){
    up <- sum(x$logFC > 0 & x$gene %in% chrX)
    down <- sum(x$logFC < 0 & x$gene %in% chrX)
    data.frame(up=up, down=down)
})
degs.chrX <- do.call(rbind, degs.chrX)
degs.chrX$down <- degs.chrX$down * -1
degs.chrX$celltype <- replace.names(gsub('_', '.', names(degs)))

degs.chrX <- melt(degs.chrX, id.vars='celltype')

# Convert variable to factor and reorder levels
degs.chrX$variable <- factor(degs.chrX$variable, levels = c("down", "up"))

pdf('Aim_1_2024/Figure_6A.pdf')
ggplot(degs.chrX, aes(x=celltype, y=value, fill=variable)) +
    geom_bar(stat="identity", position="identity") +
    coord_flip() +
    scale_y_continuous(labels = abs) +
    theme_minimal() +
    labs(x="Cell Type", y="Number of Genes", fill="Direction") +
    scale_fill_brewer(palette="Set1") +
    ggtitle('Differentially expressed chrX genes')
dev.off()

# Figure 6B - Heatmap of differentially expressed genes across cell types coloured by logFC
degs_mtx <- matrix(0, nrow=length(unique(unlist(lapply(degs, function(x) x$gene)))), ncol=length(degs))
rownames(degs_mtx) <- unique(unlist(lapply(degs, function(x) x$gene)))
colnames(degs_mtx) <- replace.names(gsub('_', '.', names(degs)))

for(i in 1:length(degs)){
    degs_mtx[degs[[i]]$gene, i] <- degs[[i]]$logFC
}

pdf('Aim_1_2024/Figure_6B.pdf')
Heatmap(t(degs_mtx), name='logFC', col=colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red')), 
        cluster_rows=TRUE, cluster_columns=TRUE, show_row_names=TRUE, show_column_names=FALSE)
dev.off()

