library(ggplot2)
library(ggrepel)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

metrics.path <- list.files('psuedobulk/ML.models/ensemble', pattern='metrics_', full.names=TRUE)

metrics <- lapply(metrics.path, function(x) {
  df <- read.csv(x)
  df <- cbind(celltype=gsub('.csv', '', basename(x)), df)
  df
})
metrics <- metrics %>% bind_rows() %>%
    mutate(perm=ifelse(grepl('perm', celltype), 'post', 'pre')) %>%
    mutate(celltype=gsub('perm_metrics_|metrics_', '', celltype)) %>%
    mutate(perm=factor(perm, levels=c('pre', 'post'))) %>%
    mutate(features=ifelse(grepl('chrX', celltype), 'chrX', 'HVG')) %>%
    mutate(features=factor(features, levels=c('chrX', 'HVG'))) %>%
    mutate(celltype=gsub('.chrX|.HVG', '', celltype)) %>%
    arrange(celltype)

# Include feature size
features.list <- list.files('psuedobulk/ML.models/ensemble/features', pattern='.txt', full.names=TRUE)
features <- lapply(features.list, function(x) {
  tmp <- read.csv(x, sep='\t', header=TRUE)
  data.frame(celltype=gsub('perm.|.chrX.txt|.HVG.txt', '', basename(x)), 
  size=nrow(tmp),
  perm=ifelse(grepl('perm.', x), 'post', 'pre'),
  features=ifelse(grepl('chrX', x), 'chrX', 'HVG'))
})
features <- bind_rows(features)

# Combine all data together
metrics$id <- paste(metrics$celltype, metrics$perm, metrics$features, sep='_')
features$id <- paste(features$celltype, features$perm, features$features, sep='_')
results <- merge(metrics, features, by='id') %>%
  select(-c(id, celltype.y, perm.y, features.y))
colnames(results) <- gsub('.x', '', colnames(results))

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')
results$celltype <- replace.names(results$celltype)

results <- rbind(results, data.frame(celltype=c('B cells', 'Age-associated B cells'), Accuracy=c(0, 0), 
Precision=c(0, 0), Recall=c(0, 0), F1=c(0, 0), F1_lower=c(0, 0), F1_upper=c(0,0), AUC=c(0,0), Kappa=c(0,0),
perm=c('post', 'post'), features=c('HVG', 'HVG'), size=c(0, 0)))

results$perm <- factor(results$perm, levels=c('pre', 'post'))
results$features <- factor(results$features, levels=c('chrX', 'HVG'))

write.table(results, 'psuedobulk/ML.models/ensemble/feature_size_metrics.txt', sep='\t', row.names=FALSE, quote=FALSE)

# boxplot of F1 score on y axis and perm on x axis split by features
pdf('psuedobulk/ML.plots/feature_perm_f1.pdf')
ggplot(metrics, aes(x=perm, y=F1, fill=features)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(0.2), alpha = 0.5) +
  facet_wrap(~features) +
  ylab('F1 score') + xlab('')
dev.off()
# Boxplot of AUC
pdf('psuedobulk/ML.plots/feature_perm_auc.pdf')
ggplot(metrics, aes(x=perm, y=AUC, fill=features)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(0.2), alpha = 0.5) +
  facet_wrap(~features) +
  ylab('AUC score') + xlab('')
dev.off()
# Boxplot of Kappa
pdf('psuedobulk/ML.plots/feature_perm_kappa.pdf')
ggplot(metrics, aes(x=perm, y=Kappa, fill=features)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(0.2), alpha = 0.5) +
  facet_wrap(~features) +
  ylab('Cohen\'s Kappa') + xlab('')
dev.off()
# Boxplot of feature size
pdf('psuedobulk/ML.plots/feature_perm_size.pdf')
ggplot(results, aes(x=perm, y=size, fill=features)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(0.2), alpha = 0.5) +
  facet_wrap(~features) +
  ylab('Feature size') + xlab('')
dev.off()


# Check normality assumptions of size to test for differences
shapiro.test(results$size[results$features=='chrX' & results$perm=='pre'])
shapiro.test(results$size[results$features=='chrX' & results$perm=='post'])
shapiro.test(results$size[results$features=='HVG' & results$perm=='pre'])
shapiro.test(results$size[results$features=='HVG' & results$perm=='post'])
# Plot histogram to view distribution
pdf('psuedobulk/ML.plots/feature_size_hist.pdf')
ggplot(results, aes(x=size, fill=perm)) +
  geom_histogram(bins=20) +
  facet_wrap(~features+perm) +
  ylab('Frequency') + xlab('Feature size')
dev.off()
### Test for differences in feature size ###
# chrX
wilcox.test(results$size[results$features=='chrX' & results$perm=='pre'], 
results$size[results$features=='chrX' & results$perm=='post'], paired=TRUE, 
alternative='greater', ties.method='FALSE', conf.int = TRUE)
# HVG
wilcox.test(results$size[results$features=='HVG' & results$perm=='pre'],
results$size[results$features=='HVG' & results$perm=='post'], paired=TRUE,
alternative='greater', ties.method='FALSE', conf.int = TRUE)
### test for differences in F1 ###
# chrX
wilcox.test(results$F1[results$features=='chrX' & results$perm=='pre'],
results$F1[results$features=='chrX' & results$perm=='post'], paired=TRUE,
alternative='greater', ties.method='FALSE', conf.int = TRUE)
# HVG
wilcox.test(results$F1[results$features=='HVG' & results$perm=='pre'],
results$F1[results$features=='HVG' & results$perm=='post'], paired=TRUE,
alternative='greater', ties.method='FALSE', conf.int = TRUE)
pdf('psuedobulk/ML.plots/feature_F1_hist.pdf')
ggplot(results, aes(x=F1)) +
  geom_histogram(bins=20) +
  facet_wrap(~features+perm) +
  ylab('Frequency') + xlab('F1 score')
dev.off()
### test for differences in AUC ###
# chrX
wilcox.test(results$AUC[results$features=='chrX' & results$perm=='pre'],
results$AUC[results$features=='chrX' & results$perm=='post'], paired=TRUE,
alternative='greater', ties.method='FALSE', conf.int = TRUE)
# HVG
wilcox.test(results$AUC[results$features=='HVG' & results$perm=='pre'],
results$AUC[results$features=='HVG' & results$perm=='post'], paired=TRUE,
alternative='greater', ties.method='FALSE', conf.int = TRUE)
### Test for differences in Kappa ###
# chrX
wilcox.test(results$Kappa[results$features=='chrX' & results$perm=='pre'],
results$Kappa[results$features=='chrX' & results$perm=='post'], paired=TRUE,
alternative='greater', ties.method='FALSE', conf.int = TRUE)
# HVG
wilcox.test(results$Kappa[results$features=='HVG' & results$perm=='pre'],
results$Kappa[results$features=='HVG' & results$perm=='post'], paired=TRUE,
alternative='greater', ties.method='FALSE', conf.int = TRUE)

# Dot plot of feature size on x axis, F1 score on y axis, fill = AUC and size = feature size
pdf('psuedobulk/ML.plots/feature_size_f1_auc.pdf')
ggplot(results, aes(y=celltype, x=F1, color=AUC)) +
    geom_point(aes(size=size)) +
    scale_color_gradient2(low = "white", mid = "blue", high = "red", 
                                             midpoint = 0.5, limits = c(0, 1), name = "AUC") +
    ylab('') + xlab('F1 score') +
    theme(axis.text.y = element_text(size=8)) +
    geom_vline(xintercept = 0.8, linetype = 'dotted') +
    facet_wrap(~features+perm) +
    theme(panel.spacing=unit(1,"lines")) +
    scale_x_continuous(limits = c(0, 1.0), expand = c(0, 0), breaks = seq(0, 1.0, by = 0.2))
dev.off()
  
# Scatter plot pre and post feature size, F1 score and AUC
results.split <- split(results, results$features)
chrX.results <- split(results.split[['chrX']], results.split[['chrX']]$perm)
chrX.results <- merge(chrX.results[['pre']], chrX.results[['post']], by='celltype', suffixes=c('.pre', '.post'))

pdf('psuedobulk/ML.plots/f1.chrX.scatter.pdf')
rho <- cor.test(chrX.results$F1.pre, chrX.results$F1.post, method='spearman', exact=FALSE)
ggplot(chrX.results, aes(x = F1.pre, y = F1.post, color = celltype)) +
  geom_point() +
  geom_text_repel(aes(label = celltype), size = 3, box.padding = 0.5, color = "black", max.overlaps = Inf) +
  theme(legend.position="none") +
  ylab('F1 score post') + xlab('F1 score pre') +
  geom_abline(intercept = 0, slope = 1, linetype='dotted') +
  xlim(0, max(chrX.results$F1.pre, chrX.results$F1.post)) +
  ylim(0, max(chrX.results$F1.pre, chrX.results$F1.post)) +
  labs(subtitle=paste('rho:', round(rho$estimate, 2), ' ', 'p-value:', 
  format(rho$p.value, scientific = TRUE, digits = 3)))
dev.off()

pdf('psuedobulk/ML.plots/auc.chrX.scatter.pdf')
rho <- cor.test(chrX.results$AUC.pre, chrX.results$AUC.post, method='spearman', exact=FALSE)
ggplot(chrX.results, aes(x = AUC.pre, y = AUC.post, color = celltype)) +
  geom_point() +
  geom_text_repel(aes(label = celltype), size = 3, box.padding = 0.5, color = "black", max.overlaps = Inf) +
  geom_abline(intercept = 0, slope = 1, linetype='dotted') +
  theme(legend.position="none") +
  ylab('AUC post') + xlab('AUC pre') +
  xlim(0, max(chrX.results$AUC.pre, chrX.results$AUC.post)) +
  ylim(0, max(chrX.results$AUC.pre, chrX.results$AUC.post)) +
  labs(subtitle=paste('rho:', round(rho$estimate, 2), ' ', 'p-value:', 
  format(rho$p.value, scientific = TRUE, digits = 3)))
dev.off()

pdf('psuedobulk/ML.plots/size.chrX.scatter.pdf')
rho <- cor.test(chrX.results$size.pre, chrX.results$size.post, method='spearman', exact=FALSE)
ggplot(chrX.results, aes(x = size.pre, y = size.post, color = celltype)) +
  geom_point() +
  geom_text_repel(aes(label = celltype), size = 3, box.padding = 0.5, color = "black", max.overlaps = Inf) +
  geom_abline(intercept = 0, slope = 1, linetype='dotted') +
  theme(legend.position="none") +
  ylab('size post') + xlab('size pre') +
  xlim(0, max(chrX.results$size.pre, chrX.results$size.post)) +
  ylim(0, max(chrX.results$size.pre, chrX.results$size.post)) +
  labs(subtitle=paste('rho:', round(rho$estimate, 2), ' ', 'p-value:', 
  format(rho$p.value, scientific = TRUE, digits = 3)))
dev.off()

HVG.results <- split(results.split[['HVG']], results.split[['HVG']]$perm)
HVG.results <- merge(HVG.results[['pre']], HVG.results[['post']], by='celltype', suffixes=c('.pre', '.post'))

pdf('psuedobulk/ML.plots/f1.HVG.scatter.pdf')
rho <- cor.test(HVG.results$F1.pre, HVG.results$F1.post, method='spearman', exact=FALSE)
ggplot(HVG.results, aes(x = F1.pre, y = F1.post, color = celltype)) +
geom_point() +
geom_text_repel(aes(label = celltype), size = 3, box.padding = 0.5, color = "black", max.overlaps = Inf) +
geom_abline(intercept = 0, slope = 1, linetype='dotted') +
xlim(0, max(HVG.results$F1.pre, HVG.results$F1.post)) +
ylim(0, max(HVG.results$F1.pre, HVG.results$F1.post)) +
theme(legend.position="none") +
ylab('F1 score post') + xlab('F1 score pre') +
labs(subtitle=paste('rho:', round(rho$estimate, 2), ' ', 'p-value:', 
  format(rho$p.value, scientific = TRUE, digits = 3)))
dev.off()

pdf('psuedobulk/ML.plots/auc.HVG.scatter.pdf')
rho <- cor.test(HVG.results$AUC.pre, HVG.results$AUC.post, method='spearman', exact=FALSE)
ggplot(HVG.results, aes(x = AUC.pre, y = AUC.post, color = celltype)) +
geom_point() +
geom_text_repel(aes(label = celltype), size = 3, box.padding = 0.8, color = "black", max.overlaps = Inf, ) +
theme(legend.position="none") +
ylab('AUC post') + xlab('AUC pre') +
geom_abline(intercept = 0, slope = 1, linetype='dotted') +
xlim(0, max(HVG.results$AUC.pre, HVG.results$AUC.post)) +
ylim(0, max(HVG.results$AUC.pre, HVG.results$AUC.post)) +
labs(subtitle=paste('rho:', round(rho$estimate, 2), ' ', 'p-value:', 
  format(rho$p.value, scientific = TRUE, digits = 3)))
dev.off()

pdf('psuedobulk/ML.plots/size.HVG.scatter.pdf')
rho <- cor.test(HVG.results$size.pre, HVG.results$size.post, method='spearman', exact=FALSE)
ggplot(HVG.results, aes(x = size.pre, y = size.post, color = celltype)) +
  geom_point() +
  geom_text_repel(aes(label = celltype), size = 3, box.padding = 0.5, color = "black", max.overlaps = Inf) +
  geom_abline(intercept = 0, slope = 1, linetype='dotted') +
  theme(legend.position="none") +
  ylab('size post') + xlab('size pre') +
  xlim(0, max(HVG.results$size.pre, HVG.results$size.post)) +
  ylim(0, max(HVG.results$size.pre, HVG.results$size.post)) +
  labs(subtitle=paste('rho:', round(rho$estimate, 2), ' ', 'p-value:', 
  format(rho$p.value, scientific = TRUE, digits = 3)))
dev.off()

# Check for X chromosome in selected features
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
perm.HVG.files <- list.files('psuedobulk/ML.models/ensemble/features', pattern='perm.*HVG', full.names=TRUE)
perm.HVG <- lapply(perm.HVG.files, read.delim)
names(perm.HVG) <- gsub('perm.|.HVG.txt', '', basename(perm.HVG.files))
lapply(perm.HVG, function(x) {
  x$Features[x$Features %in% rownames(chrX)]
})
"TSC22D3"
"RPL36A"

HVG.files <- list.files('psuedobulk/ML.models/ensemble/features', pattern='.*HVG', full.names=TRUE)
HVG <- lapply(HVG.files, read.delim)
names(HVG) <- gsub('.HVG.txt', '', basename(HVG.files))
lapply(HVG, function(x) {
  x$Features[x$Features %in% rownames(chrX)]
})

### Visualise differential expression of selected features ###
# Read in differentially expressed genes
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
degs <- deg.list('differential.expression/edgeR/', filter=FALSE)
names(degs) <- gsub('_', '.', names(degs))
deg.chrX <- unique(unlist(lapply(degs, function(x) x$gene[x$gene %in% rownames(chrX)])))

# Read in feature list - chrX pre
chrX.feature.files <- list.files('psuedobulk/ML.models/ensemble/features', pattern='chrX', full.names=TRUE)
pre.chrX.files <- chrX.feature.files[grep('perm', chrX.feature.files, invert=TRUE)]
pre.chrX <- lapply(pre.chrX.files, read.delim)
names(pre.chrX) <- gsub('.chrX.txt', '', basename(pre.chrX.files))
pre.chrX.features <- unique(unlist(lapply(pre.chrX, function(x) x$Features)))

# Create matrix
mtx <- matrix(0, nrow=length(pre.chrX.features), ncol=length(pre.chrX))
rownames(mtx) <- pre.chrX.features
colnames(mtx) <- names(pre.chrX)
# Fill matrix
for(i in 1:length(pre.chrX)){
    cell <- names(pre.chrX)[i]
    edgeR <- degs[[cell]]
    mtx[rownames(mtx) %in% edgeR$gene, i] <- ifelse(edgeR[match(rownames(mtx), edgeR$gene, nomatch=0), 'FDR'] < 0.05, 
    edgeR[match(rownames(mtx), edgeR$gene, nomatch=0), 'logFC.disease_vs_control'], 0)
}
colnames(mtx) <- replace.names(colnames(mtx))
# Plot heatmap
pdf('psuedobulk/ML.plots/pre.chrX.heatmap.pdf')
# ha = rowAnnotation(genes = anno_mark(at = match(deg.chrX, rownames(mtx)), 
#     labels = deg.chrX, labels_gp = gpar(fontsize = 10)))
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
Heatmap(mtx, column_title = "Pre chrX features", name = "logFC", col=col_fun,
clustering_distance_rows = "euclidean", clustering_distance_columns = 'euclidean',
clustering_method_rows = "complete", clustering_method_columns = "complete",
show_row_names = FALSE, column_names_gp = gpar(fontsize = 5), column_names_rot = 45)
dev.off()

# Read in feature list - chrX post
chrX.feature.files <- list.files('psuedobulk/ML.models/ensemble/features', pattern='chrX', full.names=TRUE)
post.chrX.files <- chrX.feature.files[grep('perm', chrX.feature.files)]
post.chrX <- lapply(post.chrX.files, read.delim)
names(post.chrX) <- gsub('perm.|.chrX.txt', '', basename(post.chrX.files))
post.chrX.features <- unique(unlist(lapply(post.chrX, function(x) x$Features)))

# Create matrix
mtx <- matrix(0, nrow=length(post.chrX.features), ncol=length(post.chrX))
rownames(mtx) <- post.chrX.features
colnames(mtx) <- gsub('perm.', '', names(post.chrX))
# Fill matrix
for(i in 1:length(post.chrX)){
    cell <- gsub('perm.', '', names(post.chrX)[i])
    edgeR <- degs[[cell]]
    mtx[rownames(mtx) %in% edgeR$gene, i] <- ifelse(edgeR[match(rownames(mtx), edgeR$gene, nomatch=0), 'FDR'] < 0.05, 
    edgeR[match(rownames(mtx), edgeR$gene, nomatch=0), 'logFC.disease_vs_control'], 0)
}
colnames(mtx) <- replace.names(colnames(mtx))
# Plot heatmap
pdf('psuedobulk/ML.plots/post.chrX.heatmap.pdf')
# ha = rowAnnotation(genes = anno_mark(at = match(deg.chrX, rownames(mtx)), 
#     labels = deg.chrX, labels_gp = gpar(fontsize = 10)))
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
Heatmap(mtx, column_title = "Post chrX features", name = "logFC", col=col_fun,
clustering_distance_rows = "euclidean", clustering_distance_columns = 'euclidean',
clustering_method_rows = "complete", clustering_method_columns = "complete",
show_row_names = FALSE, column_names_gp = gpar(fontsize = 5), column_names_rot = 45)
dev.off()

# Read in feature list - HVG pre
HVG.feature.files <- list.files('psuedobulk/ML.models/ensemble/features', pattern='HVG', full.names=TRUE)
pre.HVG.files <- HVG.feature.files[grep('perm', HVG.feature.files, invert=TRUE)]
pre.HVG <- lapply(pre.HVG.files, read.delim)
names(pre.HVG) <- gsub('.HVG.txt', '', basename(pre.HVG.files))
pre.HVG.features <- unique(unlist(lapply(pre.HVG, function(x) x$Features)))

# Create matrix
mtx <- matrix(0, nrow=length(pre.HVG.features), ncol=length(pre.HVG))
rownames(mtx) <- pre.HVG.features
colnames(mtx) <- names(pre.HVG)
# Fill matrix
for(i in 1:length(pre.HVG)){
    cell <- names(pre.HVG)[i]
    edgeR <- degs[[cell]]
    mtx[rownames(mtx) %in% edgeR$gene, i] <- ifelse(edgeR[match(rownames(mtx), edgeR$gene, nomatch=0), 'FDR'] < 0.05, 
    edgeR[match(rownames(mtx), edgeR$gene, nomatch=0), 'logFC.disease_vs_control'], 0)
}
colnames(mtx) <- replace.names(colnames(mtx))
# Plot heatmap
pdf('psuedobulk/ML.plots/pre.HVG.heatmap.pdf')
# ha = rowAnnotation(genes = anno_mark(at = match(deg.chrX, rownames(mtx)),
#     labels = deg.chrX, labels_gp = gpar(fontsize = 10)))
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
Heatmap(mtx, column_title = "Pre HVG features", name = "logFC", col=col_fun,
clustering_distance_rows = "euclidean", clustering_distance_columns = 'euclidean',
clustering_method_rows = "complete", clustering_method_columns = "complete",
show_row_names = FALSE, column_names_gp = gpar(fontsize = 5), column_names_rot = 45,
use_raster = TRUE)
dev.off()

# Read in feature list - HVG post
HVG.feature.files <- list.files('psuedobulk/ML.models/ensemble/features', pattern='HVG', full.names=TRUE)
post.HVG.files <- HVG.feature.files[grep('perm', HVG.feature.files)]
post.HVG <- lapply(post.HVG.files, read.delim)
names(post.HVG) <- gsub('perm.|.HVG.txt', '', basename(post.HVG.files))
post.HVG.features <- unique(unlist(lapply(post.HVG, function(x) x$Features)))

# Create matrix
mtx <- matrix(0, nrow=length(post.HVG.features), ncol=length(post.HVG))
rownames(mtx) <- post.HVG.features
colnames(mtx) <- gsub('perm.', '', names(post.HVG))
# Fill matrix
for(i in 1:length(post.HVG)){
    cell <- gsub('perm.', '', names(post.HVG)[i])
    edgeR <- degs[[cell]]
    mtx[rownames(mtx) %in% edgeR$gene, i] <- ifelse(edgeR[match(rownames(mtx), edgeR$gene, nomatch=0), 'FDR'] < 0.05, 
    edgeR[match(rownames(mtx), edgeR$gene, nomatch=0), 'logFC.disease_vs_control'], 0)
}
colnames(mtx) <- replace.names(colnames(mtx))
# Plot heatmap
pdf('psuedobulk/ML.plots/post.HVG.heatmap.pdf')
# ha = rowAnnotation(genes = anno_mark(at = match(deg.chrX, rownames(mtx)),
#     labels = deg.chrX, labels_gp = gpar(fontsize = 10)))
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
Heatmap(mtx, column_title = "Post HVG features", name = "logFC", col=col_fun,
clustering_distance_rows = "euclidean", clustering_distance_columns = 'euclidean',
clustering_method_rows = "complete", clustering_method_columns = "complete",
show_row_names = FALSE, column_names_gp = gpar(fontsize = 5), column_names_rot = 45,
use_raster = TRUE)
dev.off()

chrX.immune <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/X-linked.immune.genes.Chang.txt')
chrX.immune$X.linked.immune.genes[chrX.immune$X.linked.immune.genes %in% post.chrX$all.cells$Features]

genes <- c("AIFM1", "ATP11C", "BTK", "CXCR3", "CXorf38", "CYBB", "DDX3X", "EDA", 
"G6PD", "GAB3", "IGBP1", "IL13RA1", "IL1RAPL1", "IL2RG", "IRAK1", "KDM6A", "MSN", 
"NKAP", "NKRF", "SASH3", "SH2D1A", "TIMP1", "TLR7", "VSIG4", "WAS")

disgene <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/SLE.tsv')
disgene$Gene[disgene$Gene %in% genes]
