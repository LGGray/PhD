# Function to calculate enrichment of selected features for XCI escape genes
library(dplyr)
library(ggplot2)
library(ggsignif)
library(reshape2)
library(UpSetR)
library(rstatix)
library(ComplexHeatmap)
library(circlize)
library(stringr)
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')

# Create ML.plots directory
# dir.create('ML.plots')

# Read in ML metrics file
metrics <- read.delim('psuedobulk/metrics/metrics.combined.txt')

# Add columns for celltype, features and ML
metrics <- metrics %>% 
    mutate(celltype = gsub('.+_|.chrX|.HVG', '', model)) %>%
    mutate(features = gsub('^.*\\.', '', model)) %>%
    mutate(ML = gsub('_.+', '', model)) %>%
    arrange(celltype)
metrics$ML <- factor(metrics$ML, levels=c('logit', 'RF', 'SVM', 'GBM', 'MLP'))
metrics$features <- factor(metrics$features, levels=c('chrX', 'HVG'))

# Replace celltype names
metrics$celltype <- replace.names(metrics$celltype)

# Order metrics by F1
metrics <- metrics[order(metrics$F1, decreasing=TRUE),]

# Plot the F1 scores for each model
pdf('psuedobulk/ML.plots/F1.forest.pdf')
ggplot(metrics, aes(x=F1, y=celltype, color = ML)) +
    geom_point(size = 1, position=position_jitter(height = 0.5, seed = 42)) +
    geom_errorbarh(
        aes(xmin = F1_lower, xmax = F1_upper),
        height = 0.2,
        position=position_jitter(height = 0.5, seed = 42)) +
    geom_vline(xintercept = 0.8, linetype = 'dotted') +
    theme_bw() +
    xlab("Weighted F1-score") + ylab("") + ggtitle("Individual model performance") +
    labs(color='Features') +
    facet_wrap(~features)
dev.off()

# Plot the AUPRC scores for each model
pdf('psuedobulk/ML.plots/AUPRC.forest.pdf')
ggplot(metrics, aes(x=AUPRC, y=celltype, color = ML)) +
    geom_point(size = 1, position=position_jitter(height = 0.5, seed = 42)) +
    geom_errorbarh(
        aes(xmin = AUPRC_lower, xmax = AUPRC_upper),
        height = 0.2,
        position=position_jitter(height = 0.5, seed = 42)) +
    geom_vline(xintercept = 0.8, linetype = 'dotted') +
    theme_bw() +
    xlab("AUPRC") + ylab("") + ggtitle("Individual model performance") +
    labs(color='Features') +
    facet_wrap(~features)
dev.off()

# Compare the ranks of F1-scores between chrX and HVG for each cell type
rank.test <- lapply(unique(metrics$celltype), function(x){
    df <- subset(metrics, celltype == x)
    wilcox.test(F1 ~ features, data=df, alternative='greater')$p.value
})
names(rank.test) <- unique(metrics$celltype)

boxplot.data <- subset(metrics, celltype %in% names(rank.test)[rank.test < 0.05])
pdf('psuedobulk/ML.plots/F1.boxplot.pdf')
ggplot(boxplot.data, aes(x=features, y=F1, fill=features)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color=ML), width=0.2) +
    geom_signif(comparisons = list(c("chrX", "HVG")), test = "wilcox.test", 
    test.args=list(alternative='greater'), map_signif_level=TRUE) +
    theme_bw() +
    xlab("") + ylab("Weighted F1 score") +
    facet_wrap(~celltype)
dev.off()

# Boxplot for all celltypes
pdf('psuedobulk/ML.plots/F1.boxplot.all.pdf', width=10, height=10)
ggplot(metrics, aes(x=features, y=F1, fill=features)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color=ML), width=0.2) +
    geom_signif(comparisons = list(c("chrX", "HVG")), test = "wilcox.test", 
    test.args=list(alternative='greater'), map_signif_level=TRUE) +
    theme_bw() +
    xlab("") + ylab("Weighted F1 score") +
    facet_wrap(~celltype, ncol=6)
dev.off()

# Plot F1 scores for each ensemble model
ensemble.files <- list.files('psuedobulk/ML.models/ensemble/', pattern='metrics_.+.csv', full.names=TRUE)
ensemble.files <- grep('perm', ensemble.files, value=TRUE, invert=TRUE)
ensemble.list <- lapply(ensemble.files, read.csv)
names(ensemble.list) <- gsub('metrics_|.chrX|.HVG|.csv', '', basename(ensemble.files))
ensemble <- bind_rows(ensemble.list, .id='celltype')
ensemble$features <- str_extract(basename(ensemble.files), "chrX|HVG")
ensemble$features <- factor(ensemble$features, levels=c('chrX', 'HVG'))
# Replace names
ensemble$celltype <- replace.names(ensemble$celltype)

# Write out file
# write.table(ensemble, 'psuedobulk/ML.models/ensemble/ensemble.metrics.txt', row.names=FALSE, quote=F, sep='\t')

ensemble <- read.delim('psuedobulk/ML.models/ensemble/ensemble.metrics.txt')

pdf('psuedobulk/ML.plots/F1.forest.ensemble.pdf')
ggplot(ensemble, aes(x=F1, y=celltype, color=factor(celltype))) +
    geom_point() +
    geom_errorbarh(
        aes(xmin = F1_lower, xmax = F1_upper)) +
    geom_vline(xintercept = 0.8, linetype = 'dotted') +
    theme_bw() +
    theme(legend.position = "none",
    plot.margin = unit(c(1,1,1,0), "cm")) +
    xlab("Weighted F1-score") + ylab("") + ggtitle("Ensemble model performance") +
    facet_wrap(~features)
dev.off()

pdf('psuedobulk/ML.plots/AUPRC.forest.ensemble.pdf')
ggplot(ensemble, aes(x=AUPRC, y=celltype, color=factor(celltype))) +
    geom_point() +
    geom_errorbarh(
        aes(xmin = AUPRC_lower, xmax = AUPRC_upper)) +
    geom_vline(xintercept = 0.8, linetype = 'dotted') +
    theme_bw() +
    theme(legend.position = "none",
    plot.margin = unit(c(1,1,1,0), "cm")) +
    xlab("AUPRC score") + ylab("") + ggtitle("Ensemble model performance") +
    facet_wrap(~features)
dev.off()

pdf('psuedobulk/ML.plots/F1.ensemble.barplot.pdf')
ggplot(ensemble, aes(x=celltype, y=F1, fill=features)) +
    geom_bar(stat='identity', position='dodge') +
    geom_errorbar(
        aes(ymin = F1_lower, ymax = F1_upper),
        position=position_dodge(width=0.9), width=0.25) +
    theme_bw() +
    geom_hline(yintercept = 0.8, linetype = 'dotted') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.margin = unit(c(1,1,2,1), "cm")) +
    xlab("") + ylab("Weighted F1 score") +
    ggtitle("Ensemble model performance")
dev.off()

pdf('psuedobulk/ML.plots/AUPRC.ensemble.barplot.pdf')
ggplot(ensemble, aes(x=celltype, y=AUPRC, fill=features)) +
    geom_bar(stat='identity', position='dodge') +
    geom_errorbar(
        aes(ymin = AUPRC_lower, ymax = AUPRC_upper),
        position=position_dodge(width=0.9), width=0.25) +
    theme_bw() +
    geom_hline(yintercept = 0.8, linetype = 'dotted') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.margin = unit(c(1,1,2,1), "cm")) +
    xlab("") + ylab("AUPRC score") +
    ggtitle("Ensemble model performance")
dev.off()

# subset ensemble data for F1_lower > 0.8 and AUPRC_lower > 0.8
top_models <- subset(ensemble, F1 > 0.8 & F1_lower > 0.8 & AUPRC > 0.8 & AUPRC_lower > 0.8 & features == 'chrX')$celltype

# For top chrX models, create heatmap of selected features
feature.list <- lapply(top_models, function(x){
    read.csv(paste0('psuedobulk/features/combined_features.', gsub("/|-| ", ".", x), '.chrX.csv'))[,1]
})
names(feature.list) <- top_models
feature.mtx <- fromList(feature.list)
rownames(feature.mtx) <- unique(unlist(feature.list))

unique(unlist(lapply(feature.list, function(x){
    x[x %in% rownames(chrX)]
})))

# Plot heatmap
pdf('psuedobulk/ML.plots/feature.heatmap.top_models.chrX.pdf')
col_fun = colorRamp2(c(0, 1), c("grey", "red"))
Heatmap(as.matrix(feature.mtx), column_title = "Selected X chromosome features",
col = col_fun, clustering_distance_rows = "euclidean", clustering_distance_columns = 'euclidean',
clustering_method_rows = "complete", clustering_method_columns = "complete",
show_row_names = TRUE, column_names_gp = gpar(fontsize = 12), column_names_rot = 45,
row_names_gp = gpar(fontsize = 5), show_heatmap_legend = FALSE)
dev.off()

# Determine if top model selected features are differentially expressed
deg <- lapply(top_models, function(x){
    cell <- gsub("/|-| ", "_", x)
    edgeR <- read.delim(paste0('differential.expression/edgeR_cellCount/', cell, '.txt'))
    edgeR <- edgeR[edgeR$gene %in% feature.list[[x]],]
    return(edgeR[,c('gene', 'logFC.disease_vs_control', 'FDR')])
})
names(deg) <- replace.names(top_models)

deg.mtx <- bind_rows(deg, .id='celltype')
deg.mtx <- reshape2::dcast(deg.mtx, gene ~ celltype, value.var='logFC.disease_vs_control')
rownames(deg.mtx) <- deg.mtx$gene
deg.mtx <- deg.mtx[,-1]
deg.mtx[is.na(deg.mtx)] <- 0

# Plot heatmap of logFC
pdf('psuedobulk/ML.plots/feature.logFC.heatmap.top_models.chrX.pdf')
Heatmap(as.matrix(deg.mtx), column_title = "Selected X chromosome features",
clustering_distance_rows = "euclidean", clustering_distance_columns = 'euclidean',
clustering_method_rows = "complete", clustering_method_columns = "complete",
show_row_names = TRUE, column_names_gp = gpar(fontsize = 10), column_names_rot = 45,
row_names_gp = gpar(fontsize = 5), show_heatmap_legend = TRUE, name='logFC')
dev.off()

# Create box plot of expression levels for selected features
top_models[1]
exp <- readRDS(paste0('psuedobulk/', gsub("/|-| ", ".", x), '.chrX.RDS'))
exp <- exp[,c('class', feature.list[[x]])]
exp[2:ncol(exp)] <- scale(exp[2:ncol(exp)])
exp.wide <- reshape2::melt(exp)
exp$variable <- factor(exp$variable, levels=c('control', 'disease'))

pdf('psuedobulk/ML.plots/feature.expression.boxplot.top_models.chrX.pdf')
ggplot(exp.wide, aes(x=variable, y=value, fill=class)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    xlab("") + ylab("z-score") +
    geom_signif(comparisons = list(c('control', 'disease')), test = "wilcox.test") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(x)
dev.off()



# For top HVG models, create heatmap of selected features
top_models <- c('Tem.Effector.helper.T.cells', 'pDC', 'CD16+.NK.cells')
feature.list <- lapply(top_models, function(x){
    read.csv(paste0('psuedobulk/features/combined_features.', x, '.HVG.csv'))[,1]
})
names(feature.list) <- replace.names(top_models)
feature.mtx <- fromList(feature.list)
rownames(feature.mtx) <- unique(unlist(feature.list))

# Plot heatmap
pdf('psuedobulk/ML.plots/feature.heatmap.top_models.HVG.pdf')
col_fun = colorRamp2(c(0, 1), c("grey", "red"))
Heatmap(as.matrix(feature.mtx), column_title = "Selected HVG features",
col = col_fun, clustering_distance_rows = "euclidean", clustering_distance_columns = 'euclidean',
clustering_method_rows = "complete", clustering_method_columns = "complete",
show_row_names = TRUE, column_names_gp = gpar(fontsize = 12), column_names_rot = 45,
row_names_gp = gpar(fontsize = 5), show_heatmap_legend = FALSE)
dev.off()

# Determine if top model selected features are differentially expressed
feature.list <- lapply(top_models, function(x){
    read.csv(paste0('psuedobulk/features/combined_features.', x, '.HVG.csv'))[,1]
})
names(feature.list) <- top_models
deg <- lapply(names(feature.list), function(x){
    cell <- gsub('\\.', '_', x)
    edgeR <- read.delim(paste0('differential.expression/edgeR/', cell, '.txt'))
    edgeR <- edgeR[edgeR$gene %in% feature.list[[x]] & edgeR$FDR < 0.05,]
    return(edgeR[,c('gene', 'logFC.disease_vs_control')])
})
names(deg) <- replace.names(top_models)

deg.mtx <- bind_rows(deg, .id='celltype')
deg.mtx <- reshape2::dcast(deg.mtx, gene ~ celltype, value.var='logFC.disease_vs_control')
rownames(deg.mtx) <- deg.mtx$gene
deg.mtx <- deg.mtx[,-1]
deg.mtx[is.na(deg.mtx)] <- 0

# Plot heatmap of logFC
pdf('psuedobulk/ML.plots/feature.logFC.heatmap.top_models.HVG.pdf')
Heatmap(as.matrix(deg.mtx), column_title = "Selected HVG features",
clustering_distance_rows = "euclidean", clustering_distance_columns = 'euclidean',
clustering_method_rows = "complete", clustering_method_columns = "complete",
show_row_names = TRUE, column_names_gp = gpar(fontsize = 10), column_names_rot = 45,
row_names_gp = gpar(fontsize = 5), show_heatmap_legend = TRUE, name='logFC')
dev.off()

load('../../datasets/XCI/chrX.Rdata')
lapply(feature.list, function(x){
    x[x %in% rownames(chrX)]
})

# Calculate enrichment of XCI escape genes
# Read in feature files
top_models <- c('Tcm.Naive.helper.T.cells', 'Regulatory.T.cells', 'Non.classical.monocytes', 'Memory.B.cells')
feature.list <- lapply(top_models, function(x){
    read.csv(paste0('psuedobulk/features/combined_features.', x, '.chrX.csv'))[,1]
})
names(feature.list) <- top_models

# Calculate enrichment of XCI escape genes
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')

# Fishers test for enrichment of XCI escape genes
lapply(names(feature.list), function(x){
    gene.list <- feature.list[[x]]
    background <- readRDS(paste0('psuedobulk/', x, '.chrX.RDS'))
    background <- colnames(background)[-c(1:3)]
    background <- background[!(background %in% gene.list)]
    a <- sum(gene.list %in% rownames(escape))
    b <- sum(!(gene.list %in% rownames(escape)))
    c <- sum(background %in% rownames(escape))
    d <- sum(!(background %in% rownames(escape)))
    chisq.test(matrix(c(a, c, b, d), nrow = 2))
})


