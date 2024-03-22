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
ensemble <- subset(ensemble, celltype != "HSC/MPP")

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

pdf('psuedobulk/ML.plots/chrX_HVG_F1.dotplot.pdf')
ggplot(ensemble, aes(x=celltype, y=F1, fill=features)) +
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8)) +
    geom_errorbar(aes(ymin=F1_lower, ymax=F1_upper), width=0.2, position=position_dodge(0.8)) +  # Add error bars
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "black") +  # Add dotted line at 0.8
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('') +
    coord_flip()
dev.off()

pdf('psuedobulk/ML.plots/chrX_HVG_AUPRC.dotplot.pdf')
ggplot(ensemble, aes(x=celltype, y=AUPRC, fill=features)) +
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8)) +
    geom_errorbar(aes(ymin=AUPRC_lower, ymax=AUPRC_upper), width=0.2, position=position_dodge(0.8)) +  # Add error bars
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "black") +  # Add dotted line at 0.8
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('') +
    coord_flip()
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

### Sample size calculations ###
library(pwr)
calculate_sample_size <- function(control_group, treatment_group, power = 0.8, sig_level = 0.05) {
  # Calculate the effect size (Cohen's d)
  effect_size <- (mean(treatment_group) - mean(control_group)) / sqrt((var(control_group) + var(treatment_group)) / 2)
  
  # Calculate the sample size
  sample_size <- pwr.t.test(d = effect_size, power = power, sig.level = sig_level, type = "two.sample", alternative = "two.sided")$n
  
  return(sample_size)
}

exp <- readRDS(paste0('psuedobulk/', gsub("/|-| ", ".", top_models[4]), '.chrX.RDS'))
exp <- exp[,c('class', feature.list[[top_models[4]]])]

sample_size <- lapply(3:ncol(exp), function(i){
    control_group <- exp[exp$class == 'control', i]
    disease_group <- exp[exp$class == 'disease', i]
    calculate_sample_size(control_group, disease_group)
})
names(sample_size) <- colnames(exp)[3:ncol(exp)]
unlist(sample_size)
df <- data.frame(feature = names(sample_size), sample_size = unlist(sample_size))
# add 15% to sample size
df$updated_sample_size <- df$sample_size + (df$sample_size * 0.15)
write.table(df, 'psuedobulk/Tcm.Naive.helper.T.cells.sample_size.txt', row.names=FALSE, quote=FALSE, sep='\t')


# plot PCA1 and PCA2 for exp
pca <- prcomp(exp[,-c(1,2)], center=TRUE, scale=TRUE)
pca$x <- cbind(class=ifelse(exp$class=='control', 0, 1), pca$x)
pdf('psuedobulk/ML.plots/PCA.Tcm.Naive.helper.T.cells.pdf')
ggplot(data.frame(pca$x), aes(x=PC1, y=PC2, color=factor(class), data)) +
    geom_point() +
    theme_bw() +
    xlab("PCA1") + ylab("PCA2") +
    ggtitle(top_models[4])
dev.off()

# Plot elbow plot
pdf('psuedobulk/ML.plots/elbow.plot.pdf')
std_dev <- pca$sdev
prop_varex <- std_dev^2 / sum(std_dev^2)
# Plot the variance explained by each principal component
plot(prop_varex, xlab = "Principal Component", ylab = "Proportion of Variance Explained", type = "b")
# Add a cumulative sum line
lines(cumsum(prop_varex), type = "b", col = "red")
dev.off()

# Plot UMAP of PCA
library(umap)
umap_result <- umap(pca$x[,1:5])
colnames(umap_result$layout) <- c('UMAP1', 'UMAP2')
pdf('psuedobulk/ML.plots/UMAP.Tcm.Naive.helper.T.cells.pdf')
ggplot(data.frame(umap_result$layout), aes(x=UMAP1, y=UMAP2, color=factor(exp$class))) +
    geom_point() +
    theme_bw() +
    xlab("UMAP1") + ylab("UMAP2") +
    ggtitle(top_models[4])
dev.off()

exp[,-1] <- scale(exp[,-1])
exp_wide <- reshape2::melt(exp)
pdf('psuedobulk/ML.plots/Tcm.Naive.helper.T.cells.violin.pdf')
ggplot(exp_wide, aes(x=variable, y=value, fill=class)) +
    geom_violin() +
    theme_bw() +
    xlab("") + ylab("z-score") +
    ggtitle(top_models[4])
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

# Read in all chrX feature files
chrX.files <- list.files('psuedobulk/features/', pattern='combined_features.+chrX.csv', full.names=TRUE)
chrX.list <- lapply(chrX.files, read.csv)
names(chrX.list) <- replace.names(gsub('combined_features.|.chrX.csv', '', basename(chrX.files)))


upset.mtx <- fromList(lapply(chrX.list, function(x){x[,1]}))
rownames(upset.mtx) <- unique(unlist(lapply(chrX.list, function(x){x[,1]})))

sort(rowSums(upset.mtx[rownames(upset.mtx) %in% disgene$Gene,]))

pdf('psuedobulk/ML.plots/upset.chrX.pdf')
upset(upset.mtx, order.by='freq', nsets = length(chrX.list), nintersects = NA, show.numbers = 'no')
dev.off()

chrX.immune <- read.delim('../../datasets/XCI/X-linked.immune.genes.Chang.txt')

# Read in all HVG feature files
HVG.files <- list.files('psuedobulk/features/', pattern='combined_features.+HVG.csv', full.names=TRUE)
HVG.list <- lapply(HVG.files, read.csv)
names(HVG.list) <- replace.names(gsub('combined_features.|.HVG.csv', '', basename(HVG.files)))

upset.mtx <- fromList(lapply(HVG.list, function(x){x[,1]}))
rownames(upset.mtx) <- unique(unlist(lapply(HVG.list, function(x){x[,1]})))

sort(rowSums(upset.mtx))[rownames(chrX)]

pdf('psuedobulk/ML.plots/upset.HVG.pdf')
upset(upset.mtx, order.by='freq', nsets = length(HVG.list), nintersects = NA, show.numbers = 'no')
dev.off()