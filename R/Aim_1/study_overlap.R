library(dplyr)
library(tidyr)
library(ggplot2)
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

# Read in gene sets
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/sex_hormones.RData')
ISG <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/arazi.2019/ISG.txt')
inflammatory <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/arazi.2019/inflammatory.txt')
X.immune <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/X-linked.immune.genes.Chang.txt')

MS <- deg.list('MS_GSE193770/differential.expression/edgeR', logfc=0.1)
pSS <- deg.list('pSS_GSE157278/differential.expression/edgeR', logfc=0.1)
UC <- deg.list('UC_GSE125527/differential.expression/edgeR', logfc=0.1)
CO <- deg.list('CD_Kong/colon/differential.expression/edgeR', logfc=0.1)
names(CO)[2] <- 'CD16-_NK_cells'
TI <- deg.list('CD_Kong/TI/differential.expression/edgeR', logfc=0.1)
names(TI)[1] <- 'CD16-_NK_cells'
SLE <- deg.list('lupus_Chun/differential.expression/edgeR', logfc=0.1)
SLE_2 <- deg.list('SLE_GSE135779/differential.expression/edgeR', logfc=0.1)

# Create tables of up and downregulated genes for each cell type
MS_metrics <- bind_rows(lapply(names(MS), function(x){
    data.frame(
    study='MS',
    celltype=x,
    upregulated = sum(MS[[x]]$logFC > 0),
    upregulated_chrX = sum(MS[[x]]$logFC > 0 & MS[[x]]$gene %in% rownames(chrX)),
    downregulated = sum(MS[[x]]$logFC < 0),
    downregulated_chrX = sum(MS[[x]]$logFC < 0 & MS[[x]]$gene %in% rownames(chrX)))
}))

pSS_metrics <- bind_rows(lapply(names(pSS), function(x){
    data.frame(
    study='pSS',
    celltype=x,
    upregulated = sum(pSS[[x]]$logFC > 0),
    upregulated_chrX = sum(pSS[[x]]$logFC > 0 & pSS[[x]]$gene %in% rownames(chrX)),
    downregulated = sum(pSS[[x]]$logFC < 0),
    downregulated_chrX = sum(pSS[[x]]$logFC < 0 & pSS[[x]]$gene %in% rownames(chrX)))
}))

UC_metrics <- bind_rows(lapply(names(UC), function(x){
    data.frame(
    study='UC',
    celltype=x,
    upregulated = sum(UC[[x]]$logFC > 0),
    upregulated_chrX = sum(UC[[x]]$logFC > 0 & UC[[x]]$gene %in% rownames(chrX)),
    downregulated = sum(UC[[x]]$logFC < 0),
    downregulated_chrX = sum(UC[[x]]$logFC < 0 & UC[[x]]$gene %in% rownames(chrX)))
}))

CO_metrics <- bind_rows(lapply(names(CO), function(x){
    data.frame(
    study='CO',
    celltype=x,
    upregulated = sum(CO[[x]]$logFC > 0),
    upregulated_chrX = sum(CO[[x]]$logFC > 0 & CO[[x]]$gene %in% rownames(chrX)),
    downregulated = sum(CO[[x]]$logFC < 0),
    downregulated_chrX = sum(CO[[x]]$logFC < 0 & CO[[x]]$gene %in% rownames(chrX)))
}))

TI_metrics <- bind_rows(lapply(names(TI), function(x){
    data.frame(
    study='TI',
    celltype=x,
    upregulated = sum(TI[[x]]$logFC > 0),
    upregulated_chrX = sum(TI[[x]]$logFC > 0 & TI[[x]]$gene %in% rownames(chrX)),
    downregulated = sum(TI[[x]]$logFC < 0),
    downregulated_chrX = sum(TI[[x]]$logFC < 0 & TI[[x]]$gene %in% rownames(chrX)))
}))

SLE_metrics <- bind_rows(lapply(names(SLE), function(x){
    data.frame(
    study='SLE',
    celltype=x,
    upregulated = sum(SLE[[x]]$logFC > 0),
    upregulated_chrX = sum(SLE[[x]]$logFC > 0 & SLE[[x]]$gene %in% rownames(chrX)),
    downregulated = sum(SLE[[x]]$logFC < 0),
    downregulated_chrX = sum(SLE[[x]]$logFC < 0 & SLE[[x]]$gene %in% rownames(chrX)))
}))

# rowbind all tables
all_metrics <- bind_rows(MS_metrics, pSS_metrics, UC_metrics, CO_metrics, TI_metrics, SLE_metrics)
all_metrics$celltype <- replace.names(gsub('_', '.', all_metrics$celltype))

write.table(all_metrics, 'Aim_1/study_metrics.txt', sep='\t', row.names=FALSE, quote=FALSE)

# Plot the up and downregulated genes for each cell type coloured by study
all_metrics_long <- all_metrics %>%
  gather(key = "direction", value = "count", upregulated, downregulated) %>%
  mutate(count = ifelse(direction == "downregulated", -1 * count, count))

all_metrics_long$celltype <- factor(all_metrics_long$celltype, levels=names(celltype_colours))

all_metrics_long$logcount <- log(abs(all_metrics_long$count))
all_metrics_long$logcount <- ifelse(all_metrics_long$direction == 'downregulated', all_metrics_long$logcount * -1, all_metrics_long$logcount)
all_metrics_long$logcount[!(is.finite(all_metrics_long$logcount))] <- 0

# Plot using ggplot2
pdf('Aim_1/DEG_metrics_dotplot.pdf', width=10, height=6)
ggplot(all_metrics_long, aes(x = celltype, y = logcount, color = study)) +
  geom_point(position = position_dodge(width = 0.75), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Number of Genes (Up/Downregulated)", x = "Cell Type", title = "Differential expression by Cell Type and Study") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
  plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_color_manual(name = "Study", values = study_colours)
dev.off()

pdf('Aim_1/DEG_metrics_study_boxplot.pdf', width=10, height=6)
ggplot(all_metrics_long, aes(x=study, y=count, fill=study)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(1, 1, 1, 1, "cm")) +
    scale_fill_manual(name = "Study", values = study_colours)
dev.off()

pdf('Aim_1/DEG_metrics_celltype_boxplot.pdf', width=10, height=6)
ggplot(all_metrics_long, aes(x=celltype, y=count, fill=celltype)) +
    geom_boxplot() +
    scale_y_continuous(labels = abs) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(1, 1, 1, 1, "cm")) +
    scale_fill_manual(name = "celltype", values = celltype_colours) +
    theme(legend.position = "none")
dev.off()

# Calculate gene set enrichment in each cell type and study
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/fishers.test.degs.R')

MS_all <- deg.list('MS_GSE193770/differential.expression/edgeR', filter=FALSE)
pSS_all <- deg.list('pSS_GSE157278/differential.expression/edgeR', filter=FALSE)
UC_all <- deg.list('UC_GSE125527/differential.expression/edgeR', filter=FALSE)
CO_all <- deg.list('CD_Kong/colon/differential.expression/edgeR', filter=FALSE)
names(CO_all)[2] <- 'CD16-_NK_cells'
TI_all <- deg.list('CD_Kong/TI/differential.expression/edgeR', filter=FALSE)
names(TI_all)[1] <- 'CD16-_NK_cells'
SLE_all <- deg.list('lupus_Chun/differential.expression/edgeR', filter=FALSE)

# Enrichment of XCI escape genes
enrichment.test_escape <- list()
for(study in c('MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')){
    tmp <- lapply(get(paste0(study, '_all')), function(x) fisher.test.edgeR(x, rownames(escape), 0.1, direction='up'))
    enrichment.test_escape[[study]] <- data.frame(celltype=names(tmp), FDR=p.adjust(sapply(tmp, function(x) x$p.value), method='fdr'))
}
names(enrichment.test_escape) <- c('MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')

# Merge the tables on celltype and fill missing with NA
enrichment.test_escape_merged <- Reduce(function(x, y) merge(x, y, by='celltype', all=TRUE), enrichment.test_escape)
colnames(enrichment.test_escape_merged) <- c('celltype', 'MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')
rownames(enrichment.test_escape_merged) <- replace.names(gsub('_', '.', enrichment.test_escape_merged$celltype))
enrichment.test_escape_merged[is.na(enrichment.test_escape_merged)] <- 1
# Write the table to file
write.table(enrichment.test_escape_merged, 'Aim_1/escape_enrichment_test.txt', sep='\t', row.names=TRUE, quote=FALSE)

# Save heatmap
xcape <- Heatmap(as.matrix(enrichment.test_escape_merged[,2:ncol(enrichment.test_escape_merged)]), name='FDR', 
col=colorRamp2(c(0, 0.01, 0.05, 1), c('red', 'pink', 'blue', 'grey')), show_row_names=TRUE,
row_names_gp=gpar(fontsize=8), column_title='XCI escape')

# Enrichment of chrX genes
enrichment.test_chrX <- list()
for(study in c('MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')){
    tmp <- lapply(get(paste0(study, '_all')), function(x) fisher.test.edgeR(x, rownames(chrX), 0.1, direction='none'))
    enrichment.test_chrX[[study]] <- data.frame(celltype=names(tmp), FDR=p.adjust(sapply(tmp, function(x) x$p.value), method='fdr'))
}
names(enrichment.test_chrX) <- c('MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')

# Merge the tables on celltype and fill missing with NA
enrichment.test_chrX_merged <- Reduce(function(x, y) merge(x, y, by='celltype', all=TRUE), enrichment.test_chrX)
colnames(enrichment.test_chrX_merged) <- c('celltype', 'MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')
rownames(enrichment.test_chrX_merged) <- replace.names(gsub('_', '.', enrichment.test_chrX_merged$celltype))
enrichment.test_chrX_merged[is.na(enrichment.test_chrX_merged)] <- 1
# Write the table to file
write.table(enrichment.test_chrX_merged, 'Aim_1/chrX_enrichment_test.txt', sep='\t', row.names=TRUE, quote=FALSE)

# Save heatmap
chrx <- Heatmap(as.matrix(enrichment.test_chrX_merged[,2:ncol(enrichment.test_chrX_merged)]), name='FDR',
col=colorRamp2(c(0, 0.01, 0.05, 1), c('red', 'pink', 'blue', 'grey')), show_row_names=TRUE,
row_names_gp=gpar(fontsize=8), column_title='X chromosome')

# Enrichment of interferon stimulated genes
enrichment.test_ISG <- list()
for(study in c('MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')){
    tmp <- lapply(get(paste0(study, '_all')), function(x) fisher.test.edgeR(x, ISG$Gene, 0.1, direction='none'))
    enrichment.test_ISG[[study]] <- data.frame(celltype=names(tmp), FDR=p.adjust(sapply(tmp, function(x) x$p.value), method='fdr'))
}
names(enrichment.test_ISG) <- c('MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')

# Merge the tables on celltype and fill missing with NA
enrichment.test_ISG_merged <- Reduce(function(x, y) merge(x, y, by='celltype', all=TRUE), enrichment.test_ISG)
colnames(enrichment.test_ISG_merged) <- c('celltype', 'MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')
rownames(enrichment.test_ISG_merged) <- replace.names(gsub('_', '.', enrichment.test_ISG_merged$celltype))
enrichment.test_ISG_merged[is.na(enrichment.test_ISG_merged)] <- 1
# Write the table to file
write.table(enrichment.test_ISG_merged, 'Aim_1/ISG_enrichment_test.txt', sep='\t', row.names=TRUE, quote=FALSE)

# Save heatmap
isg <- Heatmap(as.matrix(enrichment.test_ISG_merged[,2:ncol(enrichment.test_ISG_merged)]), name='FDR',
col=colorRamp2(c(0, 0.01, 0.05, 1), c('red', 'pink', 'blue', 'grey')), show_row_names=TRUE,
row_names_gp=gpar(fontsize=8), column_title='Interferon Stimulated Genes')

# Enrichment of inflammatory genes
enrichment.test_inflammatory <- list()
for(study in c('MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')){
    tmp <- lapply(get(paste0(study, '_all')), function(x) fisher.test.edgeR(x, inflammatory$Gene, 0.1, direction='none'))
    enrichment.test_inflammatory[[study]] <- data.frame(celltype=names(tmp), FDR=p.adjust(sapply(tmp, function(x) x$p.value), method='fdr'))
}
names(enrichment.test_inflammatory) <- c('MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')

# Merge the tables on celltype and fill missing with NA
enrichment.test_inflammatory_merged <- Reduce(function(x, y) merge(x, y, by='celltype', all=TRUE), enrichment.test_inflammatory)
colnames(enrichment.test_inflammatory_merged) <- c('celltype', 'MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')
rownames(enrichment.test_inflammatory_merged) <- replace.names(gsub('_', '.', enrichment.test_inflammatory_merged$celltype))
enrichment.test_inflammatory_merged[is.na(enrichment.test_inflammatory_merged)] <- 1
# Write the table to file
write.table(enrichment.test_inflammatory_merged, 'Aim_1/inflammatory_enrichment_test.txt', sep='\t', row.names=TRUE, quote=FALSE)

# Save heatmap
inflam <- Heatmap(as.matrix(enrichment.test_inflammatory_merged[,2:ncol(enrichment.test_inflammatory_merged)]), name='FDR',
col=colorRamp2(c(0, 0.01, 0.05, 1), c('red', 'pink', 'blue', 'grey')), show_row_names=TRUE,
row_names_gp=gpar(fontsize=8), column_title='Inflammatory')

# Enrichment of AR genes
enrichment.test_AR <- list()
for(study in c('MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')){
    tmp <- lapply(get(paste0(study, '_all')), function(x) fisher.test.edgeR(x, sex_hormones$AR, 0.1, direction='none'))
    enrichment.test_AR[[study]] <- data.frame(celltype=names(tmp), FDR=p.adjust(sapply(tmp, function(x) x$p.value), method='fdr'))
}
names(enrichment.test_AR) <- c('MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')

# Merge the tables on celltype and fill missing with NA
enrichment.test_AR_merged <- Reduce(function(x, y) merge(x, y, by='celltype', all=TRUE), enrichment.test_AR)
colnames(enrichment.test_AR_merged) <- c('celltype', 'MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')
rownames(enrichment.test_AR_merged) <- replace.names(gsub('_', '.', enrichment.test_AR_merged$celltype))
enrichment.test_AR_merged[is.na(enrichment.test_AR_merged)] <- 1
# Write the table to file
write.table(enrichment.test_AR_merged, 'Aim_1/AR_enrichment_test.txt', sep='\t', row.names=TRUE, quote=FALSE)

# Save heatmap
ar <- Heatmap(as.matrix(enrichment.test_AR_merged[,2:ncol(enrichment.test_AR_merged)]), name='FDR',
col=colorRamp2(c(0, 0.01, 0.05, 1), c('red', 'pink', 'blue', 'grey')), show_row_names=TRUE,
row_names_gp=gpar(fontsize=8), column_title='Androgen Receptor')

# Enrichment of ER genes
enrichment.test_ER <- list()
for(study in c('MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')){
    tmp <- lapply(get(paste0(study, '_all')), function(x) fisher.test.edgeR(x, sex_hormones$ER, 0.1, direction='none'))
    enrichment.test_ER[[study]] <- data.frame(celltype=names(tmp), FDR=p.adjust(sapply(tmp, function(x) x$p.value), method='fdr'))
}
names(enrichment.test_ER) <- c('MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')

# Merge the tables on celltype and fill missing with NA
enrichment.test_ER_merged <- Reduce(function(x, y) merge(x, y, by='celltype', all=TRUE), enrichment.test_ER)
colnames(enrichment.test_ER_merged) <- c('celltype', 'MS', 'pSS', 'UC', 'CO', 'TI', 'SLE')
rownames(enrichment.test_ER_merged) <- replace.names(gsub('_', '.', enrichment.test_ER_merged$celltype))
enrichment.test_ER_merged[is.na(enrichment.test_ER_merged)] <- 1
# Write the table to file
write.table(enrichment.test_ER_merged, 'Aim_1/ER_enrichment_test.txt', sep='\t', row.names=TRUE, quote=FALSE)

# Save heatmap
er <- Heatmap(as.matrix(enrichment.test_ER_merged[,2:ncol(enrichment.test_ER_merged)]), name='FDR',
col=colorRamp2(c(0, 0.01, 0.05, 1), c('red', 'pink', 'blue', 'grey')), show_row_names=TRUE,
row_names_gp=gpar(fontsize=8), column_title='Estrogen Receptor')

# Plot the heatmaps
pdf('Aim_1/enrichment_heatmaps.pdf', width=20, height=12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 4)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(chrx, show_heatmap_legend = FALSE,newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(xcape, show_heatmap_legend = FALSE, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
draw(isg, show_heatmap_legend = FALSE, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(inflam, show_heatmap_legend = FALSE, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
draw(ar, show_heatmap_legend = FALSE, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
draw(er, show_heatmap_legend = FALSE, newpage = FALSE)
upViewport()
lgd = Legend(at = seq(0, 0.05, by = 0.01), col_fun = colorRamp2(c(0, 0.01, 0.05, 1), c('red', 'pink', 'blue', 'grey')), 
title = "FDR")
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
grid.draw(lgd)
upViewport()
dev.off()

### Correlation of all cell types across all studies ###
MS_merged <- bind_rows(MS, .id = 'celltype')
MS_merged$study_celltype <- paste('MS', MS_merged$celltype, sep=':')
pSS_merged <- bind_rows(pSS, .id = 'celltype')
pSS_merged$study_celltype <- paste('pSS', pSS_merged$celltype, sep=':')
UC_merged <- bind_rows(UC, .id = 'celltype')
UC_merged$study_celltype <- paste('UC', UC_merged$celltype, sep=':')
CO_merged <- bind_rows(CO, .id = 'celltype')
CO_merged$study_celltype <- paste('CO', CO_merged$celltype, sep=':')
TI_merged <- bind_rows(TI, .id = 'celltype')
TI_merged$study_celltype <- paste('TI', TI_merged$celltype, sep=':')
SLE_merged <- bind_rows(SLE, .id = 'celltype')
SLE_merged$study_celltype <- paste('SLE', SLE_merged$celltype, sep=':')

all_merged <- bind_rows(MS_merged, pSS_merged, UC_merged, CO_merged, TI_merged, SLE_merged)

# # Correlation of unfiltered data
# MS_all_merged <- bind_rows(MS_all, .id = 'celltype')
# MS_all_merged$study_celltype <- paste('MS', MS_all_merged$celltype, sep=':')
# pSS_all_merged <- bind_rows(pSS_all, .id = 'celltype')
# pSS_all_merged$study_celltype <- paste('pSS', pSS_all_merged$celltype, sep=':')
# UC_all_merged <- bind_rows(UC_all, .id = 'celltype')
# UC_all_merged$study_celltype <- paste('UC', UC_all_merged$celltype, sep=':')
# CO_all_merged <- bind_rows(CO_all, .id = 'celltype')
# CO_all_merged$study_celltype <- paste('CO', CO_all_merged$celltype, sep=':')
# TI_all_merged <- bind_rows(TI_all, .id = 'celltype')
# TI_all_merged$study_celltype <- paste('TI', TI_all_merged$celltype, sep=':')
# SLE_all_merged <- bind_rows(SLE_all, .id = 'celltype')
# SLE_all_merged$study_celltype <- paste('SLE', SLE_all_merged$celltype, sep=':')

# unfiltered_merged <- bind_rows(MS_all_merged, pSS_all_merged, UC_all_merged, CO_all_merged, TI_all_merged, SLE_all_merged)

# common_celltypes <- Reduce(intersect, list(names(pSS), names(UC), names(CO), names(TI), names(SLE)))
# all_merged <- subset(all_merged, celltype %in% common_celltypes)

### Construct matrix of all genes and celltypes ###
logFC_matrix <- all_merged[,c('gene', 'logFC', 'study_celltype')] %>% 
spread(study_celltype, logFC, fill=0)
# logFC_matrix[,-1] <- apply(logFC_matrix[,-1], 2, as.numeric)

# logFC_matrix <- subset(logFC_matrix, gene %in% rownames(escape))

# Remove genes with many NA values
pdf('Aim_1/missing_values.pdf')
hist(apply(logFC_matrix[, -1], 1, function(x) sum(is.na(x))))
dev.off()

correlation_matrix <- cor(logFC_matrix[, -1], method = "spearman", use = "pairwise.complete.obs")
correlation_matrix[is.na(correlation_matrix)] <- 0
correlation_matrix <- correlation_matrix[apply(correlation_matrix, 1, function(x) sd(x) != 0),
apply(correlation_matrix, 2, function(x) sd(x) != 0)]

# # Function to get p-value
# get_pvalue <- function(x, y) {
#   cor.test(x, y, method = "spearman")$p.value
# }
# # Apply function to each pair of columns
# pvalue_matrix <- apply(logFC_matrix[,2:ncol(logFC_matrix)], 2, function(x) {
#   apply(logFC_matrix[,2:ncol(logFC_matrix)], 2, get_pvalue, y = x)
# })
# # Replace correlation with 0 if p-value > 0.05
# correlation_matrix <- ifelse(pvalue_matrix > 0.05, 0, correlation_matrix)
# correlation_matrix[is.na(correlation_matrix)] <- 0

# Hierarchical clustering
d <- dist(correlation_matrix, method = "euclidean")
hc <- hclust(d, method = "average")

wcss <- sapply(1:10, function(k){
  kmeans_result <- kmeans(d, centers = k, nstart = 25)
  return (kmeans_result$tot.withinss)
})
pdf('Aim_1/wcss.pdf')
plot(1:10, wcss, type = "b", xlab = "Number of Clusters (k)", ylab = "Within-Cluster Sum of Squares (WCSS)")
dev.off()

sil_width <- sapply(2:10, function(k) {
  km <- kmeans(d, centers = k, nstart = 25)
  silhouette_width <- silhouette(km$cluster, d)
  mean(silhouette_width[, 3])
})
sil_width

pdf('Aim_1/correlation_matrix.pdf')
celltype = replace.names(gsub('_', '.', gsub('.+:', '', colnames(correlation_matrix))))
column_ha = HeatmapAnnotation(disease = gsub(':.+', '', rownames(correlation_matrix)), 
celltype = celltype, col = list(disease=unlist(study_colours), celltype=celltype_colours))
row_ha = rowAnnotation(celltype = celltype, disease = gsub(':.+', '', colnames(correlation_matrix)),
show_legend = c(FALSE, FALSE), col = list(celltype=celltype_colours, disease=unlist(study_colours)))
Heatmap(correlation_matrix, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), name='Rho',
show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha, right_annotation = row_ha,
clustering_distance_rows = 'euclidean', clustering_method_rows = 'average',
clustering_distance_columns = 'euclidean', clustering_method_columns = 'average',
row_km = 4, column_km = 4)
dev.off()

ht <- Heatmap(correlation_matrix, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), name='Rho',
show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha, right_annotation = row_ha,
clustering_distance_rows = 'euclidean', clustering_method_rows = 'average',
clustering_distance_columns = 'euclidean', clustering_method_columns = 'average',
row_km = 4, column_km = 4)

pdf('Aim_1/kmeans.pdf', width=10, height=6)
kmeans_clustering <- cutree(hc, k=4)
plot(hc, hang = -1, cex = 0.5)
rect.hclust(hc, k = 4, border = 2:5)
dev.off()

sil_scores <- silhouette(kmeans_clustering, d)

pdf('Aim_1/sil_scores.pdf')
plot(sil_scores)
dev.off()

### XCI escape correlation ###
logFC_matrix_escape <- all_merged %>% filter(gene %in% rownames(escape)) %>% 
select(gene, logFC, study_celltype) %>% spread(study_celltype, logFC, fill=NA)

correlation_matrix_escape <- cor(logFC_matrix_escape[, -1], method = "spearman", use = "pairwise.complete.obs")
correlation_matrix_escape[is.na(correlation_matrix_escape)] <- 0
correlation_matrix_escape <- correlation_matrix_escape[apply(correlation_matrix_escape, 1, function(x) sd(x) != 0),
apply(correlation_matrix_escape, 2, function(x) sd(x) != 0)]

pdf('Aim_1/correlation_matrix_escape.pdf')
celltype = replace.names(gsub('_', '.', gsub('.+:', '', colnames(correlation_matrix_escape))))
column_ha = HeatmapAnnotation(disease = gsub(':.+', '', rownames(correlation_matrix_escape)), 
celltype = celltype, col = list(disease=unlist(study_colours), celltype=celltype_colours))
row_ha = rowAnnotation(celltype = celltype, disease = gsub(':.+', '', colnames(correlation_matrix_escape)),
show_legend = c(FALSE, FALSE), col = list(celltype=celltype_colours, disease=unlist(study_colours)))
Heatmap(correlation_matrix_escape, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), name='Rho',
show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha, right_annotation = row_ha,
clustering_distance_rows = 'euclidean', clustering_method_rows = 'average',
clustering_distance_columns = 'euclidean', clustering_method_columns = 'average')
dev.off()

# Hierarchical clustering
d <- dist(correlation_matrix, method = "euclidean")
hc <- hclust(d, method = "average")

wcss <- sapply(1:10, function(k){
  kmeans_result <- kmeans(d, centers = k, nstart = 25)
  return (kmeans_result$tot.withinss)
})
pdf('Aim_1/wcss.pdf')
plot(1:10, wcss, type = "b", xlab = "Number of Clusters (k)", ylab = "Within-Cluster Sum of Squares (WCSS)")
dev.off()


sil_width <- sapply(2:10, function(k) {
  km <- kmeans(d, centers = k, nstart = 25)
  silhouette_width <- silhouette(km$cluster, d)
  mean(silhouette_width[, 3])
})
sil_width

pdf('Aim_1/kmeans_escape.pdf', width=10, height=6)
kmeans_clustering <- cutree(hc, k=3)
plot(hc, hang = -1, cex = 0.5)
rect.hclust(hc, k = 3, border = 2:4)
dev.off()

sil_scores <- silhouette(kmeans_clustering, d)

pdf('Aim_1/sil_scores_escape.pdf')
plot(sil_scores)
dev.off()


### Jaccard similarity ###
# Jaccard similarity
jaccard <- function(x, y) {
    x <- as.character(x)
    y <- as.character(y)
    intersect <- length(intersect(x, y))
    union <- length(union(x, y))
    return(intersect / union)
}

MS_genes <- lapply(MS, function(x) subset(x$gene))
names(MS_genes) <- paste('MS', names(MS_genes), sep=':')
pSS_genes <- lapply(pSS, function(x) x$gene)
names(pSS_genes) <- paste('pSS', names(pSS_genes), sep=':')
UC_genes <- lapply(UC, function(x) x$gene)
names(UC_genes) <- paste('UC', names(UC_genes), sep=':')
CO_genes <- lapply(CO, function(x) x$gene)
names(CO_genes) <- paste('CO', names(CO_genes), sep=':')
TI_genes <- lapply(TI, function(x) x$gene)
names(TI_genes) <- paste('TI', names(TI_genes), sep=':')
SLE_genes <- lapply(SLE, function(x) x$gene)
names(SLE_genes) <- paste('SLE', names(SLE_genes), sep=':')

# Combine all gene sets
all_genes <- c(MS_genes, pSS_genes, UC_genes, CO_genes, TI_genes, SLE_genes)

# Calculate Jaccard similarity for all combinations
results <- expand.grid(names(all_genes), names(all_genes))
colnames(results) <- c("CellType1", "CellType2")
results$JaccardIndex <- mapply(function(x, y) jaccard(all_genes[[x]], all_genes[[y]]),
                               results$CellType1, results$CellType2)
results <- spread(results, CellType2, JaccardIndex)
results[is.na(results)] <- 0

save(results, file='Aim_1/jaccard_similarity.RData')

# Replace diagonal with 1
diag(results[,-1]) <- 1

pdf('Aim_1/jaccard_similarity.heatmap.pdf')
celltype = replace.names(gsub('_', '.', gsub('.+:', '', results$CellType1)))
column_ha = HeatmapAnnotation(disease = gsub(':.+', '', colnames(results)[-1]),
celltype = celltype, col = list(disease=unlist(study_colours), celltype=celltype_colours))
row_ha = rowAnnotation(celltype = celltype, disease = gsub(':.+', '', results$CellType1),
show_legend = c(FALSE, FALSE), col = list(disease=unlist(study_colours), celltype=celltype_colours))
Heatmap(as.matrix(results[,2:ncol(results)]), col=colorRamp2(c(0, 1), c("white", "red")), name='Jaccard Index',
show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha, right_annotation = row_ha)
dev.off()

# Calculate Jaccard similarity for upregulated XCI escape genes
MS_genes <- lapply(MS, function(x) subset(x, logFC > 0 & gene %in% rownames(escape))$gene)
names(MS_genes) <- paste('MS', names(MS_genes), sep=':')
pSS_genes <- lapply(pSS, function(x) subset(x, logFC > 0 & gene %in% rownames(escape))$gene)
names(pSS_genes) <- paste('pSS', names(pSS_genes), sep=':')
UC_genes <- lapply(UC, function(x) subset(x, logFC > 0 & gene %in% rownames(escape))$gene)
names(UC_genes) <- paste('UC', names(UC_genes), sep=':')
CO_genes <- lapply(CO, function(x) subset(x, logFC > 0 & gene %in% rownames(escape))$gene)
names(CO_genes) <- paste('CO', names(CO_genes), sep=':')
TI_genes <- lapply(TI, function(x) subset(x, logFC > 0 & gene %in% rownames(escape))$gene)
names(TI_genes) <- paste('TI', names(TI_genes), sep=':')
SLE_genes <- lapply(SLE, function(x) subset(x, logFC > 0 & gene %in% rownames(escape))$gene)
names(SLE_genes) <- paste('SLE', names(SLE_genes), sep=':')

up_escape_genes <- c(MS_genes, pSS_genes, UC_genes, CO_genes, TI_genes, SLE_genes)

results_up_escape <- expand.grid(names(up_escape_genes), names(up_escape_genes))
colnames(results_up_escape) <- c("CellType1", "CellType2")
results_up_escape$JaccardIndex <- mapply(function(x, y) jaccard(up_escape_genes[[x]], up_escape_genes[[y]]),
                                      results_up_escape$CellType1, results_up_escape$CellType2)
results_up_escape <- spread(results_up_escape, CellType2, JaccardIndex)
results_up_escape[is.na(results_up_escape)] <- 0

save(results_up_escape, file='Aim_1/jaccard_similarity_up_escape.RData')

# Replace diagonal with 1
diag(results_up_escape[,-1]) <- 1

pdf('Aim_1/jaccard_similarity_up_escape.heatmap.pdf')
celltype = replace.names(gsub('_', '.', gsub('.+:', '', results_up_escape$CellType1)))
column_ha = HeatmapAnnotation(disease = gsub(':.+', '', colnames(results_up_escape)[-1]),
celltype = celltype, col = list(disease=unlist(study_colours), celltype=celltype_colours))
row_ha = rowAnnotation(celltype = celltype, disease = gsub(':.+', '', results_up_escape$CellType1),
show_legend = c(FALSE, FALSE), col = list(disease=unlist(study_colours), celltype=celltype_colours))
Heatmap(as.matrix(results_up_escape[,-1]), col=colorRamp2(c(0, 1), c("white", "red")), name='Jaccard Index',
show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha, right_annotation = row_ha,
clustering_distance_rows = 'euclidean', clustering_method_rows = 'average',
clustering_distance_columns = 'euclidean', clustering_method_columns = 'average')
dev.off()

# Hierarchical clustering
d <- dist(results_up_escape[,-1], method = "euclidean")
hc <- hclust(d, method = "average")

wcss <- sapply(1:10, function(k){
  kmeans_result <- kmeans(d, centers = k, nstart = 25)
  return (kmeans_result$tot.withinss)
})
pdf('Aim_1/wcss.pdf')
plot(1:10, wcss, type = "b", xlab = "Number of Clusters (k)", ylab = "Within-Cluster Sum of Squares (WCSS)")
dev.off()


sil_width <- sapply(2:10, function(k) {
  km <- kmeans(d, centers = k, nstart = 25)
  silhouette_width <- silhouette(km$cluster, d)
  mean(silhouette_width[, 3])
})
sil_width

pdf('Aim_1/kmeans_up_escape_.pdf', width=10, height=6)
plot(hc, hang = -1, cex = 0.5)
rect.hclust(hc, k = 2, border = 2:4)
dev.off()

kmeans_clustering <- cutree(hc, k=2)
sil_scores <- silhouette(kmeans_clustering, d)

pdf('Aim_1/sil_scores_up_escape_jaccard.pdf')
plot(sil_scores)
dev.off()

# Calculate Jaccard similarity for downregulated XCI escape genes
MS_genes <- lapply(MS, function(x) subset(x, logFC < 0 & gene %in% rownames(escape))$gene)
names(MS_genes) <- paste('MS', names(MS_genes), sep=':')
pSS_genes <- lapply(pSS, function(x) subset(x, logFC < 0 & gene %in% rownames(escape))$gene)
names(pSS_genes) <- paste('pSS', names(pSS_genes), sep=':')
UC_genes <- lapply(UC, function(x) subset(x, logFC < 0 & gene %in% rownames(escape))$gene)
names(UC_genes) <- paste('UC', names(UC_genes), sep=':')
CO_genes <- lapply(CO, function(x) subset(x, logFC < 0 & gene %in% rownames(escape))$gene)
names(CO_genes) <- paste('CO', names(CO_genes), sep=':')
TI_genes <- lapply(TI, function(x) subset(x, logFC < 0 & gene %in% rownames(escape))$gene)
names(TI_genes) <- paste('TI', names(TI_genes), sep=':')
SLE_genes <- lapply(SLE, function(x) subset(x, logFC < 0 & gene %in% rownames(escape))$gene)
names(SLE_genes) <- paste('SLE', names(SLE_genes), sep=':')

down_escape_genes <- c(MS_genes, pSS_genes, UC_genes, CO_genes, TI_genes, SLE_genes)

results_down_escape <- expand.grid(names(down_escape_genes), names(down_escape_genes))
colnames(results_down_escape) <- c("CellType1", "CellType2")
results_down_escape$JaccardIndex <- mapply(function(x, y) jaccard(down_escape_genes[[x]], down_escape_genes[[y]]),
                                      results_down_escape$CellType1, results_down_escape$CellType2)
results_down_escape <- spread(results_down_escape, CellType2, JaccardIndex)
results_down_escape[is.na(results_down_escape)] <- 0

save(results_down_escape, file='Aim_1/jaccard_similarity_down_escape.RData')

# Replace diagonal with 1
diag(results_down_escape[,-1]) <- 1

pdf('Aim_1/jaccard_similarity_down_escape.heatmap.pdf')
celltype = replace.names(gsub('_', '.', gsub('.+:', '', results_down_escape$CellType1)))
column_ha = HeatmapAnnotation(disease = gsub(':.+', '', colnames(results_down_escape)[-1]),
celltype = celltype, col = list(disease=unlist(study_colours), celltype=celltype_colours))
row_ha = rowAnnotation(celltype = celltype, disease = gsub(':.+', '', results_down_escape$CellType1),
show_legend = c(FALSE, FALSE), col = list(disease=unlist(study_colours), celltype=celltype_colours))
Heatmap(as.matrix(results_down_escape[,-1]), col=colorRamp2(c(0, 1), c("white", "red")), name='Jaccard Index',
show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha, right_annotation = row_ha,
clustering_distance_rows = 'euclidean', clustering_method_rows = 'average',
clustering_distance_columns = 'euclidean', clustering_method_columns = 'average')
dev.off()


