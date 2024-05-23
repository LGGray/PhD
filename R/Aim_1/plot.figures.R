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
X.immune <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/X.immune.txt')[,1]

### Figure 2A - UMAP coloured by cell type ###
pdf('Aim_1_2024/Figure_2A.pdf', width = 10, height = 10)
DimPlot(pbmc, group.by='cellTypist', label=FALSE, cols=celltype_colours, order=TRUE, raster=TRUE)
dev.off()

### Figure 2B - UMAP coloured by condition ###
pdf('Aim_1_2024/Figure_2B.pdf')
DimPlot(pbmc, group.by='condition', label=FALSE, cols=c('control'='#328AE3', 'disease'='#E33D32'), order=FALSE, raster=TRUE)
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
    output.logit[output.logit$PropRatio > 2.5,]
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

pdf('Aim_1_2024/Figure_5A.pdf')
ggplot(degs.df, aes(x=celltype, y=value, fill=variable)) +
    geom_bar(stat="identity", position="identity") +
    coord_flip() +
    scale_fill_manual(values=c('down'='#0B3EE6', 'up'='#E62512')) +
    # scale_y_continuous(breaks = seq(-100, 200, by = 25)) +
    theme_minimal() +
    labs(x="Cell Type", y="Number of DEGs", fill="Direction")
dev.off()

write.table(degs.df, 'Aim_1_2024/figure.data/Figure_5A.txt', sep='\t')

### Figure 5B - Representitive volcano plot ###
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

pdf('Aim_1_2024/Figure_5B.pdf')
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

# Figure 5C - UpSet plot of DEGs
deg.list.up <- fromList(lapply(degs, function(x) subset(x, logFC > 0)$gene))
colnames(deg.list.up) <- replace.names(gsub('_', '.', colnames(deg.list.up)))
pdf('Aim_1_2024/Figure_5C.pdf', onefile=F, width=10, height=10)
upset(deg.list.up, order.by = "freq", main.bar.color = "black", 
sets.bar.color = 'black', matrix.color = "black", nsets=ncol(deg.list.up),
show.numbers = 'yes')
dev.off()

# Figure 5D - UpSet plot of DEGs
deg.list.down <- fromList(lapply(degs, function(x) subset(x, logFC < 0)$gene))
colnames(deg.list.down) <- replace.names(gsub('_', '.', colnames(deg.list.down)))
pdf('Aim_1_2024/Figure_5D.pdf',onefile=F, width=10, height=10)
upset(deg.list.down, order.by = "freq", main.bar.color = "black",
sets.bar.color = 'black', matrix.color = "black", nsets=ncol(deg.list.down),
show.numbers = 'yes')
dev.off()


interaction_sizes <- sapply(deg.list.up, sum)

### Figure 6A - Barplot of up/downregulated chrX genes ###
escape <- read.delim('../../datasets/XCI/Katsir.escape.txt')
escape <- escape$Gene.Symbol
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
degs_mtx <- matrix(0, nrow=length(unique(unlist(lapply(degs.chrX, function(x) x$gene)))), ncol=length(degs.chrX))
rownames(degs_mtx) <- unique(unlist(lapply(degs.chrX, function(x) x$gene)))
colnames(degs_mtx) <- replace.names(gsub('_', '.', names(degs)))

for(i in 1:length(degs.chrX)){
    degs_mtx[degs.chrX[[i]]$gene, i] <- degs.chrX[[i]]$logFC
}

pdf('Aim_1_2024/Figure_6B.pdf')
ha = rowAnnotation(foo = anno_mark(at = which(rownames(degs_mtx) %in% X.immune), 
    labels = rownames(degs_mtx)[rownames(degs_mtx) %in% X.immune],
    labels_gp = gpar(fontsize=10)))
Heatmap(degs_mtx, name='logFC', col=colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red')), 
        cluster_rows=TRUE, cluster_columns=TRUE, 
        show_row_names=FALSE,
        show_column_names=TRUE, column_names_gp = gpar(fontsize = 12),
        right_annotation = ha)
dev.off()

### Figure 6C - Fisher's exact test of gene sets ###
edgeR <- deg.list('differential.expression/edgeR', filter=FALSE)
names(edgeR) <- gsub('CD16__NK_cells', 'CD16-_NK_cells', names(edgeR))

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
    mainPathways$escape_genes <- unlist(lapply(mainPathways$leadingEdge, function(x) sum(x %in% rownames(escape))))

    mainPathways <- data.frame(mainPathways)

    fgsea_list[[i]] <- mainPathways
}
names(fgsea_list) <- replace.names(gsub('_', '.', names(edgeR)))

save(fgsea_list, file='Aim_1_2024/figure.data/fgsea_list.RData')

load('Aim_1_2024/figure.data/fgsea_list.RData')

fgsea_df <- dplyr::bind_rows(fgsea_list, .id = "celltype")
fgsea_df$celltype <- factor(fgsea_df$celltype)

# Add a frequency column that counts each pathway occurrence
fgsea_df <- fgsea_df %>%
  group_by(pathway) %>%
  mutate(freq = n()) %>%
  mutate(celltype = forcats::fct_reorder(celltype, -log10(padj), .desc = TRUE)) %>%
  ungroup() %>%
  mutate(pathway = forcats::fct_reorder(pathway, freq, .desc = FALSE)) %>%
  arrange(desc(freq), desc(-log10(padj))) %>%
  data.frame()

pdf('lupus_Chun/Aim_1_2024/Figure_7.pdf', width = 10, height = 10)
ggplot(fgsea_df, aes(x=pathway, y=-log10(padj), fill=celltype)) + 
    geom_col(position='dodge') + 
    coord_flip() + 
    theme_minimal() +
    scale_fill_manual(values=celltype_colours) +
    xlab('Pathway') + 
    ylab('-log10(padj)') + 
    ggtitle('GSEA Hallmark pathways') + 
    theme(plot.title = element_text(size = 18, hjust = 0))
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

pdf('Aim_1_2024/Figure_10A.pdf')
ggplot(data.frame(null_distribution), aes(x=null_distribution)) + 
    geom_density(fill='blue', alpha=0.5) + 
    geom_vline(xintercept=actual_diff, linetype='dashed', color='red') + 
    theme_minimal() + 
    xlab('Difference in means') + 
    ylab('Density') + 
    ggtitle(gsub('_', ': ', auc_sig[1, 'celltype_TF']))
dev.off()

reg <- read.csv('SCENIC/reg.csv', skip=3, header=F)
reg_top_hit <- subset(reg, V1 %in% 'IRF1')

gene_targets <- lapply(reg_top_hit$V9, function(x) gsub("'", '', unlist(stringr::str_extract_all(x, "'[^']+'"))))

treg <- read.delim('differential.expression/edgeR/Regulatory_T_cells.txt')
subset(treg, gene %in% unique(unlist(gene_targets)))[,c('gene', 'logFC', 'FDR')]

plasma <- read.delim('differential.expression/edgeR/Plasma_cells.txt')
subset(plasma, gene %in% unique(unlist(gene_targets)))[,c('gene', 'logFC', 'FDR')]