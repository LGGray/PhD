library(Seurat)
pbmc <- readRDS('pbmc.female.RDS')
metadata <- pbmc@meta.data
auc_mtx <- read.csv('SCENIC/auc_mtx.csv', row.names=1)
colnames(auc_mtx) <- gsub('\\.', '', colnames(auc_mtx))
# Match cellIDs to pbmc
common_rows <- intersect(rownames(metadata), rownames(auc_mtx))
metadata <- metadata[rownames(metadata) %in% common_rows, ]
auc_mtx <- auc_mtx[rownames(auc_mtx) %in% common_rows, ]

# Add condition and celltype to auc_mtx
auc_mtx$condition <- metadata$condition[match(rownames(metadata), rownames(auc_mtx))]
auc_mtx$celltype <- metadata$cellTypist[match(rownames(metadata), rownames(auc_mtx))]

# reshape data for plotting and statistical analysis
auc_mtx_melt <- reshape2::melt(auc_mtx, id.vars = c('condition', 'celltype'))
auc_mtx_melt$condition <- factor(auc_mtx_melt$condition, levels = c('disease', 'control'))

# Plot distribution of control and disease AUC values split by TF
pdf('SCENIC/auc_distribution.pdf')
ggplot(auc_mtx_melt, aes(x = value, fill = condition)) + 
geom_density(alpha = 0.5) + 
facet_wrap(~variable, scales = 'free_y') + theme_minimal()
dev.off()

pdf('SCENIC/auc_boxplot.pdf', width = 15, height = 15)
ggplot(auc_mtx_melt, aes(x = variable, y = value, fill = condition)) + 
    geom_boxplot(outlier.shape = NA) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~celltype, scales = 'free_y')   
dev.off()

auc_mtx_melt$celltype.TF <- paste(auc_mtx_melt$celltype, auc_mtx_melt$variable, sep = ':')

# Wilcoxon rank sum test
auc_mtx_split <- split(auc_mtx_melt, auc_mtx_melt$celltype.TF)
wilcox.two_sided <- lapply(auc_mtx_split, function(x) {
    wilcox.test(value ~ condition, data = x)$p.value
})

wilcox.greater <- lapply(auc_mtx_split, function(x) {
  wilcox.test(value ~ condition, data = x, alternative = 'greater')$p.value
})

wilcox.less <- lapply(auc_mtx_split, function(x) {
  wilcox.test(value ~ condition, data = x, alternative = 'less')$p.value
})

wilcox <- data.frame(two_sided = unlist(wilcox.two_sided), greater = unlist(wilcox.greater), less = unlist(wilcox.less))
# Split celltype.TF into celltype and TF columns
wilcox$celltype <- gsub(':.*', '', rownames(wilcox))
wilcox$TF <- gsub('.*:', '', rownames(wilcox))

# Bonferroni correction
wilcox$two_sided <- p.adjust(wilcox$two_sided, method = 'fdr')
wilcox$greater <- p.adjust(wilcox$greater, method = 'fdr')
wilcox$less <- p.adjust(wilcox$less, method = 'fdr')

# Save results
write.csv(wilcox, 'SCENIC/wilcox_test.csv')

library(ComplexHeatmap)
library(circlize)
library(ggplot2)

pdf('SCENIC/wilcox_heatmap.pdf')
row_ha = rowAnnotation(TF=wilcox$TF, celltype = wilcox$celltype, show_legend = c(TRUE, TRUE))
Heatmap(as.matrix(wilcox[,1:3]), name = 'p-value', col = colorRamp2(c(0.05, 0.5, 1), c('red', 'white', 'blue')), 
show_row_names = FALSE, right_annotation = row_ha)
dev.off()

# Greater in disease
wilcox_greater <- wilcox[wilcox$greater < 0.05,]