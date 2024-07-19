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
# load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')
escape <- read.delim('../../datasets/XCI/Katsir.escape.txt')
escape <- escape$Gene.Symbol
disgene <- read.delim(paste0('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/', disease, '.tsv'))$Gene
GWAS <- read.delim(paste0('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/', disease, '_GWAS.tsv'))
GWAS <- unique(unlist(lapply(GWAS$Gene, function(x) unlist(strsplit(x, ';')))))

chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name

chrY <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrY_biomaRt.txt')
chrY <- subset(chrY, Gene.name != '')
chrY <- chrY$Gene.name
# /load('/directflow/SCCGGroupShare/projects/lacgra/datasets/sex_hormones.RData')
X.immune <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/X.immune.txt')[,1]

hallmark <- clusterProfiler::read.gmt('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/h.all.v7.5.1.symbols.gmt')
positional <- clusterProfiler::read.gmt('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/c1.all.v2023.1.Hs.symbols.gmt')
estrogen <- unique(subset(hallmark, term %in% c('HALLMARK_ESTROGEN_RESPONSE_EARLY','HALLMARK_ESTROGEN_RESPONSE_LATE'))$gene)
androgen <- unique(subset(hallmark, term %in% c('HALLMARK_ANDROGEN_RESPONSE'))$gene)

### Figure 1 - UMAP coloured by cell type ###
pdf('Aim_1_2024/Figure_1.pdf', width = 10, height = 10)
DimPlot(pbmc, group.by='cellTypist', label=FALSE, cols=celltype_colours, order=TRUE, raster=TRUE)
dev.off()

### Figure 2 - UMAP coloured by condition ###
pdf('Aim_1_2024/Figure_2.pdf')
DimPlot(pbmc, group.by='condition', label=FALSE, cols=c('control'='#328AE3', 'disease'='#E33D32'), order=FALSE, raster=TRUE)
dev.off()

### Figure S1 - Marker gene expression dotplot ###
load('/directflow/SCCGGroupShare/projects/lacgra/PhD/R/marker_genes.RData')
average_expression <- AverageExpression(pbmc, features = marker_genes)$decontXcounts
marker_genes <- rownames(average_expression[rowSums(average_expression) > 0,])

pdf('Aim_1_2024/Figure_S1.pdf', width = 20, height = 10)
DotPlot(pbmc, features = marker_genes, cols = c("lightgrey", "red")) + 
theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8)) +
ylab('') + xlab('Marker genes')
dev.off()

### Figure S2 - Cell type proportions clustered bar chart ###
props <- getTransformedProps(clusters = pbmc$cellTypist, 
                             sample = pbmc$individual)

pdf('Aim_1_2024/Figure_S2.pdf', width = 10, height = 10)
p <- plotCellTypeProps(clusters = pbmc$cellTypist, sample = pbmc$individual) + 
    ggtitle("Cell type proportions") + 
    theme(plot.title = element_text(size = 18, hjust = 0)) +
    scale_fill_manual(values = celltype_colours)
p
dev.off()

write.table(p$data, 'Aim_1_2024/figure.data/Figure_S2.txt', sep='\t')

### Figure 4B - Cell type proportions scatter plot ###
pbmc$condition <- factor(pbmc$condition, levels=c('disease', 'control'))
output.logit <- propeller(clusters=pbmc$cellTypist, sample=pbmc$individual, 
    group=pbmc$condition, transform='logit')

# Create a new data frame based on the condition
label_data <- if(nrow(output.logit[output.logit$FDR < 0.05,]) > 0){
    output.logit[output.logit$FDR < 0.05,]
} else {
    output.logit[output.logit$PropRatio > 2.5,]
}

pdf('Aim_1_2024/Figure_3.pdf')
ggplot(output.logit, aes(x=log2(PropRatio), y=-log10(FDR), colour=BaselineProp.clusters)) + 
    geom_point() +
    scale_colour_manual(values=celltype_colours) +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    geom_vline(xintercept=log2(1), linetype='dashed') +
    geom_text_repel(data=label_data, aes(label=BaselineProp.clusters), 
                    box.padding = 0.5, point.padding = 0.5, colour='black') +
    theme_minimal() + 
    theme(legend.position = 'none') + 
    xlab('Proportion Ratio') + 
    ylab('-log10(FDR)') + 
    ggtitle('Cell type proportions') + 
    theme(plot.title = element_text(size = 18, hjust = 0))
dev.off()

write.table(output.logit, 'Aim_1_2024/figure.data/Figure_3.txt', sep='\t')
write.table(label_data, 'Aim_1_2024/figure.data/Figure_3_label.txt', sep='\t')

# Plot the proportions of each cell type
output.logit <- output.logit[order(output.logit$BaselineProp.Freq, decreasing=FALSE),]
output.logit$BaselineProp.clusters <- factor(output.logit$BaselineProp.clusters, levels = output.logit$BaselineProp.clusters)
pdf('Aim_1_2024/Celltype_prop.pdf')
ggplot(output.logit, aes(x=BaselineProp.clusters, y=BaselineProp.Freq, fill=BaselineProp.clusters)) + 
    geom_bar(stat='identity') +
    scale_fill_manual(values=celltype_colours) + 
    coord_flip() + 
    theme_minimal() + 
    theme(legend.position = 'none') + 
    xlab('Cell type') + 
    ylab('Proportion') + 
    ggtitle('Cell type proportions') + 
    theme(plot.title = element_text(size = 18, hjust = 0))
dev.off()


### Figure 4 - Barplot of up/downregulated genes ###
degs <- deg.list('differential.expression/edgeR', logfc=0.1)
names(degs) <- gsub('CD16__NK_cells', 'CD16-_NK_cells', names(degs))
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

degs.df <- degs.df %>%
  group_by(celltype) %>%
  mutate(total_degs = sum(abs(value))) %>%
  ungroup() %>%
  arrange(total_degs)
# Convert celltype to factor
degs.df$celltype <- factor(degs.df$celltype, levels = unique(degs.df$celltype))

pdf('Aim_1_2024/Figure_4.pdf')
ggplot(degs.df, aes(x=celltype, y=value, fill=variable)) +
    geom_bar(stat="identity", position="identity") +
    coord_flip() +
    scale_fill_manual(values=c('down'='#0B3EE6', 'up'='#E62512')) +
    theme(axis.text.y = element_text(size=18)) +
    theme_minimal() +
    labs(x="", y="Number of genes", fill="Direction", title="Differentially expressed genes")
dev.off()

write.table(degs.df, 'Aim_1_2024/figure.data/Figure_4.txt', sep='\t')

### Figure 5 - Representitive volcano plot ###
top_celltype <- names(sort(unlist(lapply(degs, nrow)), decreasing=TRUE)[1])
example <- read.delim(paste0('differential.expression/edgeR/', top_celltype, '.txt'), sep='\t')

example$colour <- ifelse(example$FDR < 0.05 & example$logFC > 0.1, "#E62512",
                            ifelse(example$FDR < 0.05 & example$logFC < -0.1, "#0B3EE6", "black"))

# Select top 10 upregulated and downregulated genes
top_up <- example %>% 
    filter(FDR < 0.05 & logFC > 0.1) %>%
    arrange(desc(logFC)) %>%
    head(n=10)
top_down <- example %>% 
    filter(FDR < 0.05 & logFC < -0.1) %>%
    arrange(logFC) %>%
    head(n=10)
top_genes <- rbind(top_up, top_down)

pdf('Aim_1_2024/Figure_5.pdf')
ggplot(example, aes(x=logFC, y=-log10(FDR), color=colour)) + 
    geom_point(alpha=0.5) +
    scale_color_identity() + 
    geom_hline(yintercept=-log10(0.05), linetype='dashed') + 
    geom_vline(xintercept=c(-0.1, 0.1), linetype='dashed') + 
    geom_text_repel(data = top_genes, aes(label = gene, color='black'), 
    size = 3, max.overlaps = Inf, nudge_x = 0.5, nudge_y = 0.5) +
    theme_minimal() + 
    theme(legend.position = 'none') + 
    xlab('logFC') + 
    ylab('-log10(FDR)') + 
    ggtitle(replace.names(gsub('_', '.', top_celltype))) + 
    theme(plot.title = element_text(size = 18, hjust = 0))
dev.off()

# Figure X - UpSet plot of all degs
deg.list.all <- fromList(lapply(degs, function(x) x$gene))
rownames(deg.list.all) <- unique(unlist(lapply(degs, function(x) x$gene)))
colnames(deg.list.all) <- replace.names(gsub('_', '.', colnames(deg.list.all)))
pdf('Aim_1_2024/UpSet_all.pdf', onefile=F, width=10, height=10)
upset(deg.list.all, order.by = "freq", main.bar.color = "black",
sets.bar.color = 'black', matrix.color = "black", nsets=ncol(deg.list.all),
show.numbers = 'yes')
dev.off()

# Figure 6 - UpSet plot of upregulated DEGs
deg.list.up <- fromList(lapply(degs, function(x) subset(x, logFC > 0)$gene))
rownames(deg.list.up) <- unique(unlist(lapply(degs, function(x) subset(x, logFC > 0)$gene)))
colnames(deg.list.up) <- replace.names(gsub('_', '.', colnames(deg.list.up)))
pdf('Aim_1_2024/Figure_6.pdf', onefile=F, width=10, height=10)
upset(deg.list.up, order.by = "freq", main.bar.color = "black", 
sets.bar.color = 'black', matrix.color = "black", nsets=ncol(deg.list.up),
show.numbers = 'yes')
dev.off()

# Figure 7 - scatter plot of Classical and Non-classical monocytes upregulated DEGs
selected_genes <- rownames(deg.list.up[rowSums(deg.list.up[, -which(colnames(deg.list.up) %in% c("Classical_monocytes", "Non_classical_monocytes"))]) == 0 
& deg.list.up$Classical_monocytes == 1 & deg.list.up$Non_classical_monocytes == 1, ])

classical <- subset(degs[['Classical_monocytes']], logFC > 0 & gene %in% selected_genes)
nonclassical <- subset(degs[['Non_classical_monocytes']], logFC > 0 & gene %in% selected_genes)
merged <- merge(classical, nonclassical, by='gene')
correlation <- cor.test(merged$logFC.x, merged$logFC.y, method='spearman')
merged$rank_logFC.x <- rank(merged$logFC.x)
merged$rank_logFC.y <- rank(merged$logFC.y)
merged <- merged[order(merged$rank_logFC.x, decreasing=TRUE),]
top_genes <- merged[1:10,]

pdf('Aim_1_2024/Figure_7.pdf')
ggplot(merged, aes(x=rank_logFC.x, y=rank_logFC.y)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, color = "red") +  
    #geom_text_repel(data = top_genes, aes(label = gene), size = 3) +
    theme_minimal() + 
    xlab('Classical monocytes ranked logFC') + 
    ylab('Non-classical monocytes ranked logFC') +
    labs(title = "Classical Monocytes vs Non-classical monocytes", 
    subtitle=(paste('Rho:', round(correlation$estimate, 2)))) +
    theme(plot.title = element_text(size = 18, hjust = 0))
dev.off()

# Figure 8 - UpSet plot of DEGs
deg.list.down <- fromList(lapply(degs, function(x) subset(x, logFC < 0)$gene))
rownames(deg.list.down) <- unique(unlist(lapply(degs, function(x) subset(x, logFC < 0)$gene)))
colnames(deg.list.down) <- replace.names(gsub('_', '.', colnames(deg.list.down)))
pdf('Aim_1_2024/Figure_8.pdf',onefile=F, width=10, height=10)
upset(deg.list.down, order.by = "freq", main.bar.color = "black",
sets.bar.color = 'black', matrix.color = "black", nsets=ncol(deg.list.down),
show.numbers = 'yes')
dev.off()

# Figure 9 - Scater plot of MAIT cells and pDC
selected_genes <- rownames(deg.list.down[rowSums(deg.list.down[, -which(colnames(deg.list.down) %in% c("MAIT_cells", "pDC"))]) == 0 
& deg.list.down$MAIT_cells == 1 & deg.list.down$pDC == 1, ])
mait <- subset(degs[['MAIT_cells']], logFC < 0 & gene %in% selected_genes)
pdc <- subset(degs[['pDC']], logFC < 0 & gene %in% selected_genes)
merged <- merge(mait, pdc, by='gene')
correlation <- cor.test(merged$logFC.x, merged$logFC.y, method='spearman')
merged$rank_logFC.x <- rank(merged$logFC.x)
merged$rank_logFC.y <- rank(merged$logFC.y)
merged <- merged[order(merged$rank_logFC.x, decreasing=FALSE),]
top_genes <- merged[1:10,]

pdf('Aim_1_2024/Figure_9.pdf')
ggplot(merged, aes(x=rank_logFC.x, y=rank_logFC.y)) + 
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, color = "red") +  
    #geom_text_repel(data = top_genes, aes(label = gene), size = 3) +
    theme_minimal() + 
    xlab('MAIT cells ranked logFC') + 
    ylab('pDC ranked logFC') +
    labs(title = "MAIT cells vs pDC", 
    subtitle=(paste('Rho:', round(correlation$estimate, 2)))) +
    theme(plot.title = element_text(size = 18, hjust = 0))
dev.off()

### Figure 6A - Barplot of up/downregulated chrX genes ###
degs.chrX <- lapply(degs, function(x){
    up.chrX.genes <- subset(x, logFC > 0 & gene %in% chrX)$gene
    up.chrX <- sum(!up.chrX.genes %in% escape)
    up.escape <- sum(up.chrX.genes %in% escape)

    down.chrX.genes <- subset(x, logFC < 0 & gene %in% chrX)$gene
    down.chrX <- sum(!down.chrX.genes %in% escape)
    down.escape <- sum(down.chrX.genes %in% escape)

    data.frame(up.chrX=up.chrX, up.escape=up.escape, down.chrX=down.chrX, down.escape=down.escape)
})
degs.chrX <- do.call(rbind, degs.chrX)
degs.chrX$down.chrX <- degs.chrX$down.chrX * -1
degs.chrX$down.escape <- degs.chrX$down.escape * -1
degs.chrX$celltype <- replace.names(gsub('_', '.', names(degs)))

degs.chrX <- melt(degs.chrX, id.vars='celltype')

# Convert variable to factor and reorder levels
degs.chrX$variable <- factor(degs.chrX$variable, levels = c("down.escape", "down.chrX", "up.escape", "up.chrX"))

degs.chrX <- degs.chrX %>%
  group_by(celltype) %>%
  mutate(total_degs = sum(abs(value))) %>%
  ungroup() %>%
  arrange(total_degs)
# Convert celltype to factor
degs.chrX$celltype <- factor(degs.chrX$celltype, levels = unique(degs.chrX$celltype))

pdf('Aim_1_2024/Figure_6A.pdf')
ggplot(degs.chrX, aes(x=celltype, y=value, fill=variable)) +
    geom_bar(stat="identity", position="stack") +
    coord_flip() +
    theme_minimal() +
    labs(x="Cell Type", y="Number of Genes", fill="Direction") +
    scale_fill_manual(values=c("down.chrX"="#0B3EE6", "down.escape"="#3EB2E5", "up.chrX"="#E62512", "up.escape"="#8F0C1E")) +
    ggtitle('Differentially expressed chrX genes')
dev.off()

### Figure 6B - Heatmap of differentially expressed genes across cell types coloured by logFC ###
degs.chrX <- lapply(degs, function(x){subset(x, gene %in% chrX)})
degs_mtx <- matrix(0, nrow=length(unique(unlist(lapply(degs.chrX, function(x) x$gene)))), ncol=length(degs.chrX))
rownames(degs_mtx) <- unique(unlist(lapply(degs.chrX, function(x) x$gene)))
colnames(degs_mtx) <- replace.names(gsub('_', '.', names(degs)))

for(i in 1:length(degs.chrX)){
    degs_mtx[degs.chrX[[i]]$gene, i] <- degs.chrX[[i]]$logFC
}

pdf('Aim_1_2024/Figure_6B.pdf')
ha = rowAnnotation(foo = anno_mark(at = which(rownames(degs_mtx) %in% rowlabels), 
    labels = rownames(degs_mtx)[rownames(degs_mtx) %in% X.immune],
    labels_gp = gpar(fontsize=10)))
Heatmap(degs_mtx, name='logFC', col=colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red')), 
        cluster_rows=TRUE, cluster_columns=TRUE, 
        show_row_names=FALSE,
        show_column_names=TRUE, column_names_gp = gpar(fontsize = 8),
        right_annotation = ha)
dev.off()

### Figure 6C - Fisher's exact test of gene sets ###
edgeR <- deg.list('differential.expression/edgeR', filter=FALSE)
names(edgeR) <- gsub('CD16__NK_cells', 'CD16-_NK_cells', names(edgeR))

chrY <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrY_biomaRt.txt')
chrY <- subset(chrY, Gene.name != '')
chrY <- chrY$Gene.name

# chrY
chrY.enrichment <- lapply(edgeR, function(x){
    tmp <- fisher.test.edgeR(x, logfc = 0.1, chrY)
    tmp$p.value
})


# XCI escape
xcape.enrichment <- lapply(edgeR, function(x){
    tmp <- fisher.test.edgeR(x, logfc = 0.1, escape)
    tmp$p.value
})
names(xcape.enrichment) <- replace.names(gsub('_', '.', names(xcape.enrichment)))
xcape.enrichment <- data.frame(
    celltype=names(xcape.enrichment), 
    p.value=unlist(xcape.enrichment), 
    FDR=p.adjust(unlist(xcape.enrichment), method='fdr')),
    size=unlist(lapply(edgeR, function(x) nrow(subset(x, FDR < 0.05 & abs(logFC) > 0.1 & gene %in% escape))))

xcape.enrichment <- xcape.enrichment[order(xcape.enrichment$size, decreasing=TRUE),]
xcape.enrichment$celltype <- factor(xcape.enrichment$celltype, levels = xcape.enrichment$celltype)
pdf('Aim_1_2024/Figure_6C.pdf')
ggplot(xcape.enrichment, aes(x=celltype, y=-log10(FDR), size=size, colour=-log10(FDR))) + 
    geom_point() +
    scale_color_gradient(low='blue', high='red') +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    xlab('') + 
    ylab('-log10(FDR)') + 
    ggtitle('Enrichment of XCI escape genes') +
    coord_flip() +
    theme(plot.margin = margin(1, 1, 3, 1))
dev.off()

# Figure X - ranked logFC values coloured by escape
cd16_nk <- edgeR[['CD16+_NK_cells']]
cd16_nk$escape <- ifelse(cd16_nk$gene %in% escape & abs(cd16_nk$logFC) > 0.1, 'Escape', 'Non-escape')
cd16_nk$rank_logFC <- rank(cd16_nk$logFC)
cd16_nk <- cd16_nk[order(cd16_nk$rank_logFC),]
pdf('Aim_1_2024/CD16+_NK_cells_ranked_logFC.pdf')
ggplot(cd16_nk, aes(x=rank_logFC, y=logFC)) +
    geom_point(aes(color = ifelse(abs(logFC) > 0.1, '#228b22', 'black')), alpha=0.3) +
    geom_text_repel(data=subset(cd16_nk, escape == 'Escape'), 
                    aes(label=gene), color='black', size=3, max.overlaps = Inf) +
    theme_minimal() +
    geom_hline(yintercept=0, linetype='dashed') + 
    xlab('Ranked logFC') + 
    ylab('logFC') + 
    ggtitle('CD16+ NK cells') + 
    theme(plot.title = element_text(size = 18, hjust = 0))
dev.off()

# Upset plot of escape genes
degs.escape <- lapply(degs, function(x){subset(x, gene %in% escape)})
deg.list.escape <- fromList(lapply(degs.chrX, function(x) x$gene))
rownames(deg.list.chrX) <- unique(unlist(lapply(degs.chrX, function(x) x$gene)))
colnames(deg.list.chrX) <- replace.names(gsub('_', '.', colnames(deg.list.chrX)))
pdf('Aim_1_2024/Figure_6D.pdf', onefile=F, width=10, height=10)
upset(deg.list.chrX, order.by = "freq", main.bar.color = "black",
sets.bar.color = 'black', matrix.color = "black", nsets=ncol(deg.list.chrX),
show.numbers = 'yes')
dev.off()

# chrX 
chrX.enrichment <- lapply(edgeR, function(x){
    tmp <- fisher.test.edgeR(x, logfc = 0.1, chrX)
    tmp$p.value
})
names(chrX.enrichment) <- replace.names(gsub('_', '.', names(chrX.enrichment)))
chrX.enrichment <- data.frame(
    celltype=names(chrX.enrichment), 
    p.value=unlist(chrX.enrichment), 
    FDR=p.adjust(unlist(chrX.enrichment), method='fdr'))

# DisGeNet
disgene.enrichment <- lapply(edgeR, function(x){
    tmp <- fisher.test.edgeR(x, logfc = 0.1, disgene)
    tmp$p.value
})
names(disgene.enrichment) <- replace.names(gsub('_', '.', names(disgene.enrichment)))
disgene.enrichment <- data.frame(
    celltype=names(disgene.enrichment), 
    p.value=unlist(disgene.enrichment),
    FDR=p.adjust(unlist(disgene.enrichment), method='fdr'),
    size=unlist(lapply(edgeR, function(x) nrow(subset(x, FDR < 0.05 & abs(logFC) > 0.1 & gene %in% disgene))))
)

# GWAS
GWAS.enrichment <- lapply(edgeR, function(x){
    tmp <- fisher.test.edgeR(x, logfc = 0.1, GWAS)
    tmp$p.value
})
names(GWAS.enrichment) <- replace.names(gsub('_', '.', names(GWAS.enrichment)))
GWAS.enrichment <- data.frame(
    celltype=names(GWAS.enrichment), 
    p.value=unlist(GWAS.enrichment),
    FDR=p.adjust(unlist(GWAS.enrichment), method='fdr'),
    size=unlist(lapply(edgeR, function(x) nrow(subset(x, FDR < 0.05 & abs(logFC) > 0.1 & gene %in% GWAS))))
)

# Estrogen
ER.enrichment <- lapply(edgeR, function(x){
    tmp <- fisher.test.edgeR(x, logfc = 0.1, estrogen)
    tmp$p.value
})
names(ER.enrichment) <- replace.names(gsub('_', '.', names(ER.enrichment)))
ER.enrichment <- data.frame(
    celltype=names(ER.enrichment), 
    p.value=unlist(ER.enrichment),
    FDR=p.adjust(unlist(ER.enrichment), method='fdr'),
    size=unlist(lapply(edgeR, function(x) nrow(subset(x, FDR < 0.05 & abs(logFC) > 0.1 & gene %in% estrogen))))
)

# Androgen
AR.enrichment <- lapply(edgeR, function(x){
    tmp <- fisher.test.edgeR(x, logfc = 0.1, androgen)
    tmp$p.value
})
names(AR.enrichment) <- replace.names(gsub('_', '.', names(AR.enrichment)))
AR.enrichment <- data.frame(
    celltype=names(AR.enrichment), 
    p.value=unlist(AR.enrichment),
    FDR=p.adjust(unlist(AR.enrichment), method='fdr'),
    size=unlist(lapply(edgeR, function(x) nrow(subset(x, FDR < 0.05 & abs(logFC) > 0.1 & gene %in% androgen))))
)

# Plot enrichment for each gene set
disgene.enrichment <- disgene.enrichment[order(disgene.enrichment$size, decreasing=TRUE),]
disgene.enrichment$celltype <- factor(disgene.enrichment$celltype, levels = disgene.enrichment$celltype)
pdf('Aim_1_2024/Disgene.enrichment.pdf')
ggplot(disgene.enrichment, aes(x=celltype, y=-log10(FDR), size=size, colour=-log10(FDR))) + 
    geom_point() +
    scale_color_gradient(low='blue', high='red') +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    xlab('') + 
    ylab('-log10(FDR)') + 
    ggtitle('Enrichment of SLE genes') +
    coord_flip() +
    theme(plot.margin = margin(1, 1, 3, 1))
dev.off()

GWAS.enrichment <- GWAS.enrichment[order(GWAS.enrichment$size, decreasing=TRUE),]
GWAS.enrichment$celltype <- factor(GWAS.enrichment$celltype, levels = GWAS.enrichment$celltype)
pdf('Aim_1_2024/GWAS.enrichment.pdf')
ggplot(GWAS.enrichment, aes(x=celltype, y=-log10(FDR), size=size, colour=-log10(FDR))) + 
    geom_point() +
    scale_color_gradient(low='blue', high='red') +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    xlab('') + 
    ylab('-log10(FDR)') + 
    ggtitle('Enrichment of GWAS genes') +
    coord_flip() +
    theme(plot.margin = margin(1, 1, 3, 1))
dev.off()

ER.enrichment <- ER.enrichment[order(ER.enrichment$size, decreasing=TRUE),]
ER.enrichment$celltype <- factor(ER.enrichment$celltype, levels = ER.enrichment$celltype)
pdf('Aim_1_2024/ER.enrichment.pdf')
ggplot(ER.enrichment, aes(x=celltype, y=-log10(FDR), size=size, colour=-log10(FDR))) + 
    geom_point() +
    scale_color_gradient(low='blue', high='red') +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    xlab('') + 
    ylab('-log10(FDR)') + 
    ggtitle('Enrichment of estrogen genes') +
    coord_flip() +
    theme(plot.margin = margin(1, 1, 3, 1))
dev.off()

AR.enrichment <- AR.enrichment[order(AR.enrichment$size, decreasing=TRUE),]
AR.enrichment$celltype <- factor(AR.enrichment$celltype, levels = AR.enrichment$celltype)
pdf('Aim_1_2024/AR.enrichment.pdf')
ggplot(AR.enrichment, aes(x=celltype, y=-log10(FDR), size=size, colour=-log10(FDR))) + 
    geom_point() +
    scale_color_gradient(low='blue', high='red') +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    xlab('') + 
    ylab('-log10(FDR)') + 
    ggtitle('Enrichment of androgen genes') +
    coord_flip() +
    theme(plot.margin = margin(1, 1, 3, 1))
dev.off()

# Biallelically expressed genes
biallelic <- lapply(degs, function(x){subset(x, gene %in% chrX & logFC < -1)$gene})
sort(table(unlist(biallelic)))
biallelic.df <- bind_rows(biallelic, .id='celltype')
biallelic_list <- fromList(biallelic)
rownames(biallelic_list) <- unique(unlist(biallelic))
pdf('Aim_1_2024/Biallelic_genes_upset.pdf', onefile=F, width=10, height=10)
upset(data.frame(t(biallelic_list)), order.by = "freq", main.bar.color = "black",
sets.bar.color = 'black', matrix.color = "black", nsets=ncol(biallelic_list),
show.numbers = 'yes')
dev.off()

### GSEA hallmark ###
hallmark <- gmtPathways('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/h.all.v7.5.1.symbols.gmt')
position <- gmtPathways('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/c1.all.v2023.1.Hs.symbols.gmt')

hallmark_fgsea_list <- list()
for(i in 1:length(edgeR)){
    # edgeR[[i]] <- subset(edgeR[[i]], !(gene %in% chrY))
    ranked_genes <- edgeR[[i]]$logFC
    names(ranked_genes) <- edgeR[[i]]$gene

    fgseaRes <- fgsea(pathways = hallmark, 
                    stats    = ranked_genes,
                    minSize  = 15,
                    maxSize  = 500)

    collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.05], 
                                        hallmark, ranked_genes)
    mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
                            order(-NES),]
    # mainPathways$escape_genes <- unlist(lapply(mainPathways$leadingEdge, function(x) sum(x %in% escape)))

    mainPathways <- data.frame(mainPathways)

    hallmark_fgsea_list[[i]] <- mainPathways
}
names(hallmark_fgsea_list) <- replace.names(gsub('_', '.', names(edgeR)))

save(hallmark_fgsea_list, file='Aim_1_2024/figure.data/fgsea_list.RData')

position_fgsea_list <- list()
for(i in 1:length(edgeR)){
    # edgeR[[i]] <- subset(edgeR[[i]], !(gene %in% chrY))
    ranked_genes <- edgeR[[i]]$logFC
    names(ranked_genes) <- edgeR[[i]]$gene

    fgseaRes <- fgsea(pathways = position, 
                    stats    = ranked_genes,
                    minSize  = 15,
                    maxSize  = 500)

    collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.05], 
                                        position, ranked_genes)
    mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
                            order(-NES),]
    # mainPathways$escape_genes <- unlist(lapply(mainPathways$leadingEdge, function(x) sum(x %in% escape)))

    mainPathways <- data.frame(mainPathways)

    position_fgsea_list[[i]] <- mainPathways
}
names(position_fgsea_list) <- replace.names(gsub('_', '.', names(edgeR)))

save(position_fgsea_list, file='Aim_1_2024/figure.data/fgsea_position.RData')

fgsea_df <- dplyr::bind_rows(fgsea_list, .id = "celltype")
fgsea_df$celltype <- factor(fgsea_df$celltype)
fgsea_df$pathway <- gsub('HALLMARK_', '', fgsea_df$pathway)

# Add a frequency column that counts each pathway occurrence
fgsea_df <- fgsea_df %>%
  group_by(pathway) %>%
  mutate(freq = n()) %>%
  ungroup() %>%
  data.frame()

# Order the pathways by frequency
fgsea_df <- fgsea_df[order(fgsea_df$freq, decreasing = FALSE),]
fgsea_df$pathway <- factor(fgsea_df$pathway, levels = unique(fgsea_df$pathway))

pdf('Aim_1_2024/Figure_7.pdf', width=10, height=10)
ggplot(fgsea_df, aes(x=celltype, y=pathway, color=NES)) + 
    geom_point() + 
    theme_minimal() +
    scale_color_gradient(low = "blue", high = "red") +
    xlab('') + 
    ylab('') + 
    ggtitle('MSigDB Hallmark pathways') + 
    theme(plot.title = element_text(size = 18, hjust = 0),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=14),
          axis.text.y = element_text(size = 14))
dev.off()

fgsea_df$escape_genes <- unlist(lapply(fgsea_df$leadingEdge, function(x) sum(x %in% X.immune)))
fgsea_df <- fgsea_df[order(fgsea_df$escape_genes, decreasing=TRUE),]

lapply(subset(fgsea_df, pathway == 'HALLMARK_ALLOGRAFT_REJECTION')[,'leadingEdge'], function(x) x[x %in% X.immune])

### KEGG Gene ontology
kegg <- gmtPathways('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/c2.cp.kegg.v7.5.1.symbols.gmt')
kegg_list <- list()
for(i in 1:length(edgeR)){
    ranked_genes <- edgeR[[i]]$logFC
    names(ranked_genes) <- edgeR[[i]]$gene

    fgseaRes <- fgsea(pathways = kegg, 
                    stats    = ranked_genes,
                    minSize  = 15,
                    maxSize  = 500)

    collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.05], 
                                        kegg, ranked_genes)
    mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
                            order(-NES),]
    mainPathways$escape_genes <- unlist(lapply(mainPathways$leadingEdge, function(x) sum(x %in% escape)))

    mainPathways <- data.frame(mainPathways)

    kegg_list[[i]] <- mainPathways
}
names(kegg_list) <- replace.names(gsub('_', '.', names(edgeR)))

save(kegg_list, file='Aim_1_2024/figure.data/go_list.RData')

kegg_df <- dplyr::bind_rows(kegg_list, .id = "celltype")
kegg_df$celltype <- factor(kegg_df$celltype)
kegg_df$pathway <- gsub('KEGG_', '', kegg_df$pathway)

# Add a frequency column that counts each pathway occurrence
kegg_df <- kegg_df %>%
  group_by(pathway) %>%
  mutate(freq = n()) %>%
  ungroup() %>%
  data.frame()

# Order the pathways by frequency
kegg_df <- kegg_df[order(kegg_df$freq, decreasing = FALSE),]
kegg_df$pathway <- factor(kegg_df$pathway, levels = unique(kegg_df$pathway))

pdf('Aim_1_2024/KEGG_GSEA.pdf')
ggplot(kegg_df, aes(x=celltype, y=pathway, color=NES)) + 
    geom_point() + 
    theme_minimal() +
    scale_color_gradient(low = "blue", high = "red") +
    xlab('') + 
    ylab('') + 
    ggtitle('GSEA KEGG ') + 
    theme(plot.title = element_text(size = 18, hjust = 0),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
dev.off()


### Figure 9 - Heatmap of SCENIC TF AUC comparison FDR ###
auc <- read.csv('SCENIC/permutation_test_avg.csv')
auc$celltype <- gsub('_.+', '', auc$celltype)
auc$TF <- gsub('.+_', '', auc$TF)

auc_mtx <- reshape2::dcast(auc, TF ~ celltype, value.var='FDR', drop=FALSE, fill=1)
rownames(auc_mtx) <- auc_mtx$TF

pdf('Aim_1_2024/Figure_9.pdf')
Heatmap(as.matrix(auc_mtx[,-1]), name='FDR', col=colorRamp2(c(0, 0.05, 1), c('red', 'blue', 'white')), 
        cluster_rows=TRUE, cluster_columns=TRUE, show_row_names=TRUE, row_names_gp = gpar(fontsize = 12), 
        show_column_names=TRUE, column_names_gp = gpar(fontsize = 12))
dev.off()

pdf('Aim_1_2024/Figure_9.pdf')
ha = rowAnnotation(foo = anno_mark(at = which(rownames(auc_mtx[,-1]) %in% chrX), 
    labels = rownames(auc_mtx[,-1])[rownames(auc_mtx[,-1]) %in% chrX],
    labels_gp = gpar(fontsize=10)))
Heatmap(as.matrix(auc_mtx[,-1]), name='FDR', col=colorRamp2(c(0, 0.05, 1), c('red', 'blue', 'white')), 
        cluster_rows=TRUE, cluster_columns=TRUE, show_row_names=FALSE, 
        show_column_names=TRUE, column_names_gp = gpar(fontsize = 12), right_annotation = ha)
dev.off()

auc_sig <- subset(auc, FDR < 0.05)
auc_sig <- auc_sig[order(auc_sig$actual_diff, decreasing=TRUE),]
auc_sig_chrX <- subset(auc, FDR < 0.05 & TF %in% chrX)
auc_sig_chrX <- auc_sig_chrX[order(auc_sig_chrX$actual_diff, decreasing=TRUE),]

### Figure 10 - distribution of significant TFs ###
metadata <- pbmc@meta.data
auc_mtx <- read.csv('SCENIC/SCENIC.auc.csv', row.names=1)
colnames(auc_mtx) <- gsub('\\.', '', colnames(auc_mtx))
# Match cellIDs to pbmc
common_rows <- intersect(rownames(metadata), rownames(auc_mtx))
metadata <- metadata[rownames(metadata) %in% common_rows, ]
auc_mtx <- auc_mtx[rownames(auc_mtx) %in% common_rows, ]

# Add condition and celltype to auc_mtx
auc_mtx$condition <- metadata$condition[match(rownames(metadata), rownames(auc_mtx))]
auc_mtx$celltype <- metadata$cellTypist[match(rownames(metadata), rownames(auc_mtx))]
auc_mtx$individual <- metadata$individual[match(rownames(metadata), rownames(auc_mtx))]

# reshape data for plotting and statistical analysis
auc_mtx_melt <- reshape2::melt(auc_mtx, id.vars = c('condition', 'celltype', 'individual'))
auc_mtx_melt$condition <- factor(auc_mtx_melt$condition, levels = c('disease', 'control'))

auc_avg <- aggregate(auc_mtx_melt$value, by = list(auc_mtx_melt$condition, auc_mtx_melt$individual, auc_mtx_melt$celltype, auc_mtx_melt$variable), FUN = mean)
colnames(auc_avg) <- c("condition", "individual", "celltype", "TF", "AUC")
auc_avg$celltype_TF <- paste(auc_avg$celltype, auc_avg$TF, sep = '_')

# Split by celltype_TF
auc_avg_split <- split(auc_avg, auc_avg$celltype_TF)

top_TF <- auc_avg_split[[auc_sig[1, 'celltype_TF']]]
# top_TF <- auc_avg_split[['Age-associated B cells_FOS']]

set.seed(123)
control_values <- top_TF[top_TF$condition == 'control','AUC']
disease_values <- top_TF[top_TF$condition == 'disease', 'AUC']

# Calculate actual difference in means
actual_diff <- mean(disease_values) - mean(control_values)

# Initialize a null distribution
null_distribution <- numeric(10000)

# Perform permutations
for (i in 1:10000) {
    # Shuffle the AUC values
    shuffled_values <- sample(c(control_values, disease_values))
    
    # Split into new groups
    perm_control <- shuffled_values[1:length(control_values)]
    perm_disease <- shuffled_values[(length(control_values)+1):length(shuffled_values)]
    
    # Calculate the difference in means for the permuted groups
    null_distribution[i] <- mean(perm_disease) - mean(perm_control)
}

# Calculate the p-value
p_value <- sum(null_distribution >= actual_diff) / length(null_distribution)

pdf('Aim_1_2024/ABC_FOS.pdf')
ggplot(data.frame(null_distribution), aes(x=null_distribution)) + 
    geom_density(fill='blue', alpha=0.5) + 
    geom_vline(xintercept=actual_diff, linetype='dashed', color='red') + 
    theme_minimal() + 
    xlab('Difference in means') + 
    ylab('Density') + 
    ggtitle('Age-associated B cells: FOS')
dev.off()

reg <- read.csv('SCENIC/reg.csv', skip=3, header=F)

hits <- lapply(sig_TF, function(TF){
    reg_top_hit <- subset(reg, V1 %in% TF)
    gene_targets <- lapply(reg_top_hit$V9, function(x) gsub("'", '', unlist(stringr::str_extract_all(x, "'[^']+'"))))
    foo <- unique(unlist(gene_targets))
    foo[foo %in% chrX]
})
names(hits) <- sig_TF

TF_fromlist <- fromList(hits)
rownames(TF_fromlist) <- unique(unlist(hits))

pdf('Aim_1_2024/TF_targets.pdf')
Heatmap(as.matrix(TF_fromlist), name='TF', col=colorRamp2(c(0, 1), c('white', 'red')), 
        cluster_rows=TRUE, cluster_columns=TRUE, show_row_names=TRUE, row_names_gp = gpar(fontsize = 12), 
        show_column_names=TRUE, column_names_gp = gpar(fontsize = 12))
dev.off()

TF_targets <- rownames(TF_fromlist)

TF_targets_deg <- lapply(degs, function(x) subset(x, gene %in% TF_targets)[,c('gene', 'logFC', 'FDR')])
plot.data <- dplyr::bind_rows(TF_targets_deg, .id = "celltype")
plot.data$celltype <- replace.names(gsub('_', '.', plot.data$celltype))
plot.data <- reshape2::dcast(plot.data, gene ~ celltype, value.var='logFC', fill=0)
rownames(plot.data) <- plot.data$gene
pdf('Aim_1_2024/TF_targets_deg.pdf')
Heatmap(as.matrix(plot.data[,-1]), name='logFC', col=colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red')), 
        cluster_rows=TRUE, cluster_columns=TRUE, show_row_names=TRUE, row_names_gp = gpar(fontsize = 12), 
        show_column_names=TRUE, column_names_gp = gpar(fontsize = 12))
dev.off()

reg_top_hit <- subset(reg, V1 %in% 'BCL11A')

gene_targets <- lapply(reg_top_hit$V9, function(x) gsub("'", '', unlist(stringr::str_extract_all(x, "'[^']+'"))))
foo <- unique(unlist(gene_targets))
hits <- foo[foo %in% chrX]

TF_targets <- lapply(degs, function(x) subset(x, gene %in% hits)[,c('gene', 'logFC', 'FDR')])
plot.data <- dplyr::bind_rows(TF_targets, .id = "celltype")
plot.data$celltype <- replace.names(gsub('_', '.', plot.data$celltype))
plot.data$celltype <- factor(plot.data$celltype)
plot.data$gene <- factor(plot.data$gene)

pdf('Aim_1_2024/CEBPB_targets.pdf')
ggplot(plot.data, aes(x=celltype, y=logFC)) +
geom_boxplot() +
geom_jitter(aes(color=gene), width=0.2) +
geom_hline(yintercept=0, linetype='dashed') +
xlab('') +
coord_flip()
dev.off()

treg <- read.delim('differential.expression/edgeR/Regulatory_T_cells.txt')
subset(treg, gene %in% unique(unlist(gene_targets)))[,c('gene', 'logFC', 'FDR')]

plasma <- read.delim('differential.expression/edgeR/Plasma_cells.txt')
subset(plasma, gene %in% unique(unlist(gene_targets)))[,c('gene', 'logFC', 'FDR')]