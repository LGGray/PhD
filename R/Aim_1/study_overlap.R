library(dplyr)
library(tidyr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')

# Jaccard similarity
jaccard <- function(x, y) {
    x <- as.character(x)
    y <- as.character(y)
    intersect <- length(intersect(x, y))
    union <- length(union(x, y))
    return(intersect / union)
}

# Read in gene sets
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')
ISG <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/arazi.2019/ISG.txt')
inflammatory <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/arazi.2019/inflammatory.txt')

MS <- deg.list('MS_GSE193770/differential.expression/edgeR', logfc=0.1)
pSS <- deg.list('pSS_GSE157278/differential.expression/edgeR', logfc=0.1)
UC <- deg.list('UC_GSE125527/differential.expression/edgeR', logfc=0.1)
CD_colon <- deg.list('CD_Kong/colon/differential.expression/edgeR', logfc=0.1)
names(CD_colon)[2] <- 'CD16-_NK_cells'
CD_TI <- deg.list('CD_Kong/TI/differential.expression/edgeR', logfc=0.1)
names(CD_TI)[1] <- 'CD16-_NK_cells'
SLE <- deg.list('lupus_Chun/differential.expression/edgeR', logfc=0.1)

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

CD_colon_metrics <- bind_rows(lapply(names(CD_colon), function(x){
    data.frame(
    study='CD_colon',
    celltype=x,
    upregulated = sum(CD_colon[[x]]$logFC > 0),
    upregulated_chrX = sum(CD_colon[[x]]$logFC > 0 & CD_colon[[x]]$gene %in% rownames(chrX)),
    downregulated = sum(CD_colon[[x]]$logFC < 0),
    downregulated_chrX = sum(CD_colon[[x]]$logFC < 0 & CD_colon[[x]]$gene %in% rownames(chrX)))
}))

CD_TI_metrics <- bind_rows(lapply(names(CD_TI), function(x){
    data.frame(
    study='CD_TI',
    celltype=x,
    upregulated = sum(CD_TI[[x]]$logFC > 0),
    upregulated_chrX = sum(CD_TI[[x]]$logFC > 0 & CD_TI[[x]]$gene %in% rownames(chrX)),
    downregulated = sum(CD_TI[[x]]$logFC < 0),
    downregulated_chrX = sum(CD_TI[[x]]$logFC < 0 & CD_TI[[x]]$gene %in% rownames(chrX)))
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
all_metrics <- bind_rows(MS_metrics, pSS_metrics, UC_metrics, CD_colon_metrics, CD_TI_metrics, SLE_metrics)
all_metrics$celltype <- replace.names(gsub('_', '.', all_metrics$celltype))

write.table(all_metrics, 'Aim_1/study_metrics.txt', sep='\t', row.names=FALSE, quote=FALSE)

# Calculate enrichment of XCI escape genes and chrX genes in each cell type and study
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/fishers.test.degs.R')

MS_all <- deg.list('MS_GSE193770/differential.expression/edgeR', filter=FALSE)
pSS_all <- deg.list('pSS_GSE157278/differential.expression/edgeR', filter=FALSE)
UC_all <- deg.list('UC_GSE125527/differential.expression/edgeR', filter=FALSE)
CD_colon_all <- deg.list('CD_Kong/colon/differential.expression/edgeR', filter=FALSE)
names(CD_colon_all)[2] <- 'CD16-_NK_cells'
CD_TI_all <- deg.list('CD_Kong/TI/differential.expression/edgeR', filter=FALSE)
names(CD_TI_all)[1] <- 'CD16-_NK_cells'
SLE_all <- deg.list('lupus_Chun/differential.expression/edgeR', filter=FALSE)

# Enrichment of XCI escape genes
enrichment.test_escape <- list()
for(study in c('MS', 'pSS', 'UC', 'CD_colon', 'CD_TI', 'SLE')){
    tmp <- lapply(get(paste0(study, '_all')), function(x) fisher.test.edgeR(x, rownames(escape), 0.1, direction='none'))
    enrichment.test_escape[[study]] <- data.frame(celltype=names(tmp), pvalue=sapply(tmp, function(x) x$p.value))
}
names(enrichment.test_escape) <- c('MS', 'pSS', 'UC', 'CD_colon', 'CD_TI', 'SLE')

# Merge the tables on celltype and fill missing with NA
enrichment.test_escape_merged <- Reduce(function(x, y) merge(x, y, by='celltype', all=TRUE), enrichment.test_escape)
colnames(enrichment.test_escape_merged) <- c('celltype', 'MS', 'pSS', 'UC', 'CD_colon', 'CD_TI', 'SLE')
rownames(enrichment.test_escape_merged) <- replace.names(gsub('_', '.', enrichment.test_escape_merged$celltype))
enrichment.test_escape_merged[is.na(enrichment.test_escape_merged)] <- 1
# Write the table to file
write.table(enrichment.test_escape_merged, 'Aim_1/escape_enrichment_test.txt', sep='\t', row.names=TRUE, quote=FALSE)

# Save heatmap
xcape <- Heatmap(as.matrix(enrichment.test_escape_merged[,2:ncol(enrichment.test_escape_merged)]), name='p-value', 
col=colorRamp2(c(0, 0.05, 1), c('red', 'blue', 'grey')), show_row_names=TRUE,
row_names_gp=gpar(fontsize=8), column_title='XCI escape')

# Enrichment of chrX genes
enrichment.test_chrX <- list()
for(study in c('MS', 'pSS', 'UC', 'CD_colon', 'CD_TI', 'SLE')){
    tmp <- lapply(get(paste0(study, '_all')), function(x) fisher.test.edgeR(x, rownames(chrX), 0.1, direction='none'))
    enrichment.test_chrX[[study]] <- data.frame(celltype=names(tmp), pvalue=sapply(tmp, function(x) x$p.value))
}
names(enrichment.test_chrX) <- c('MS', 'pSS', 'UC', 'CD_colon', 'CD_TI', 'SLE')

# Merge the tables on celltype and fill missing with NA
enrichment.test_chrX_merged <- Reduce(function(x, y) merge(x, y, by='celltype', all=TRUE), enrichment.test_chrX)
colnames(enrichment.test_chrX_merged) <- c('celltype', 'MS', 'pSS', 'UC', 'CD_colon', 'CD_TI', 'SLE')
rownames(enrichment.test_chrX_merged) <- replace.names(gsub('_', '.', enrichment.test_chrX_merged$celltype))
enrichment.test_chrX_merged[is.na(enrichment.test_chrX_merged)] <- 1
# Write the table to file
write.table(enrichment.test_chrX_merged, 'Aim_1/chrX_enrichment_test.txt', sep='\t', row.names=TRUE, quote=FALSE)

# Save heatmap
chrx <- Heatmap(as.matrix(enrichment.test_chrX_merged[,2:ncol(enrichment.test_chrX_merged)]), name='p-value',
col=colorRamp2(c(0, 0.05, 1), c('red', 'blue', 'grey')), show_row_names=TRUE,
row_names_gp=gpar(fontsize=8), column_title='X chromosome')

# Enrichment of interferon stimulated genes
enrichment.test_ISG <- list()
for(study in c('MS', 'pSS', 'UC', 'CD_colon', 'CD_TI', 'SLE')){
    tmp <- lapply(get(paste0(study, '_all')), function(x) fisher.test.edgeR(x, ISG$Gene, 0.1, direction='none'))
    enrichment.test_ISG[[study]] <- data.frame(celltype=names(tmp), pvalue=sapply(tmp, function(x) x$p.value))
}
names(enrichment.test_ISG) <- c('MS', 'pSS', 'UC', 'CD_colon', 'CD_TI', 'SLE')

# Merge the tables on celltype and fill missing with NA
enrichment.test_ISG_merged <- Reduce(function(x, y) merge(x, y, by='celltype', all=TRUE), enrichment.test_ISG)
colnames(enrichment.test_ISG_merged) <- c('celltype', 'MS', 'pSS', 'UC', 'CD_colon', 'CD_TI', 'SLE')
rownames(enrichment.test_ISG_merged) <- replace.names(gsub('_', '.', enrichment.test_ISG_merged$celltype))
enrichment.test_ISG_merged[is.na(enrichment.test_ISG_merged)] <- 1
# Write the table to file
write.table(enrichment.test_ISG_merged, 'Aim_1/ISG_enrichment_test.txt', sep='\t', row.names=TRUE, quote=FALSE)

# Save heatmap
isg <- Heatmap(as.matrix(enrichment.test_ISG_merged[,2:ncol(enrichment.test_ISG_merged)]), name='p-value',
col=colorRamp2(c(0, 0.05, 1), c('red', 'blue', 'grey')), show_row_names=TRUE,
row_names_gp=gpar(fontsize=8), column_title='ISG')

# Enrichment of inflammatory genes
enrichment.test_inflammatory <- list()
for(study in c('MS', 'pSS', 'UC', 'CD_colon', 'CD_TI', 'SLE')){
    tmp <- lapply(get(paste0(study, '_all')), function(x) fisher.test.edgeR(x, inflammatory$Gene, 0.1, direction='none'))
    enrichment.test_inflammatory[[study]] <- data.frame(celltype=names(tmp), pvalue=sapply(tmp, function(x) x$p.value))
}
names(enrichment.test_inflammatory) <- c('MS', 'pSS', 'UC', 'CD_colon', 'CD_TI', 'SLE')

# Merge the tables on celltype and fill missing with NA
enrichment.test_inflammatory_merged <- Reduce(function(x, y) merge(x, y, by='celltype', all=TRUE), enrichment.test_inflammatory)
colnames(enrichment.test_inflammatory_merged) <- c('celltype', 'MS', 'pSS', 'UC', 'CD_colon', 'CD_TI', 'SLE')
rownames(enrichment.test_inflammatory_merged) <- replace.names(gsub('_', '.', enrichment.test_inflammatory_merged$celltype))
enrichment.test_inflammatory_merged[is.na(enrichment.test_inflammatory_merged)] <- 1
# Write the table to file
write.table(enrichment.test_inflammatory_merged, 'Aim_1/inflammatory_enrichment_test.txt', sep='\t', row.names=TRUE, quote=FALSE)

# Save heatmap
inflam <- Heatmap(as.matrix(enrichment.test_inflammatory_merged[,2:ncol(enrichment.test_inflammatory_merged)]), name='p-value',
col=colorRamp2(c(0, 0.05, 1), c('red', 'blue', 'grey')), show_row_names=TRUE,
row_names_gp=gpar(fontsize=8), column_title='Inflammatory')

# Plot the heatmaps
pdf('Aim_1/enrichment_heatmaps.pdf', width=12, height=12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 3)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(xcape, show_heatmap_legend = FALSE,newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(chrx, show_heatmap_legend = FALSE, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(isg, show_heatmap_legend = FALSE, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
draw(inflam, show_heatmap_legend = FALSE, newpage = FALSE)
upViewport()
lgd = Legend(at = c(-2, -1, 0, 1, 2), col_fun = colorRamp2(c(0, 0.05, 1), c('red', 'blue', 'grey')), 
title = "pvalue")
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
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
CD_colon_merged <- bind_rows(CD_colon, .id = 'celltype')
CD_colon_merged$study_celltype <- paste('CD_colon', CD_colon_merged$celltype, sep=':')
CD_TI_merged <- bind_rows(CD_TI, .id = 'celltype')
CD_TI_merged$study_celltype <- paste('CD_TI', CD_TI_merged$celltype, sep=':')
SLE_merged <- bind_rows(SLE, .id = 'celltype')
SLE_merged$study_celltype <- paste('SLE', SLE_merged$celltype, sep=':')

all_merged <- bind_rows(MS_merged, pSS_merged, UC_merged, CD_colon_merged, CD_TI_merged, SLE_merged)

common_celltypes <- Reduce(intersect, list(names(pSS), names(UC), names(CD_colon), names(CD_TI), names(SLE)))
all_merged <- subset(all_merged, celltype %in% common_celltypes)

library(RColorBrewer)
cols_study = brewer.pal(6, "Spectral")
names(cols_study) = c('CD_colon', 'CD_TI', 'MS', 'pSS', 'SLE', 'UC')

load('/directflow/SCCGGroupShare/projects/lacgra/PhD/R/celltype.colours.RData')
cols_celltype <- colours[replace.names(gsub('_', '.', common_celltypes))]

### All genes ###
logFC_matrix <- all_merged[,c('gene', 'logFC', 'study_celltype')] %>% 
spread(study_celltype, logFC, fill=0)
correlation_matrix <- cor(logFC_matrix[,2:ncol(logFC_matrix)], method = "spearman")

pdf('Aim_1/correlation_matrix.pdf')
celltype = replace.names(gsub('_', '.', gsub('.+:', '', colnames(correlation_matrix))))
column_ha = HeatmapAnnotation(disease = gsub(':.+', '', rownames(correlation_matrix)), 
celltype = celltype, col = list(disease=cols_study, celltype=cols_celltype))
row_ha = rowAnnotation(disease = gsub(':.+', '', colnames(correlation_matrix)),
celltype = celltype, show_legend = c(FALSE, FALSE), col = list(disease=cols_study, celltype=cols_celltype))
Heatmap(correlation_matrix, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), name='Rho',
show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha, right_annotation = row_ha)
dev.off()

### chrX correlation ###
logFC_matrix_chrX <- all_merged %>% filter(gene %in% rownames(chrX)) %>% select(gene, logFC, study_celltype) %>% spread(study_celltype, logFC, fill=0)
correlation_matrix_chrX <- cor(logFC_matrix_chrX[,2:ncol(logFC_matrix_chrX)], method = "spearman")

pdf('Aim_1/correlation_matrix_chrX.pdf')
celltype = replace.names(gsub('_', '.', gsub('.+:', '', colnames(correlation_matrix_chrX))))
column_ha = HeatmapAnnotation(disease = gsub(':.+', '', rownames(correlation_matrix_chrX)),
celltype = celltype, col = list(disease=cols_study, celltype=cols_celltype))
row_ha = rowAnnotation(disease = gsub(':.+', '', colnames(correlation_matrix_chrX)),
celltype = celltype, show_legend = c(FALSE, FALSE), col = list(disease=cols_study, celltype=cols_celltype))
Heatmap(correlation_matrix_chrX, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), name='Rho',
show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha, right_annotation = row_ha)
dev.off()

### XCI escape correlation ###
logFC_matrix_escape <- all_merged %>% filter(gene %in% rownames(escape)) %>% 
select(gene, logFC, study_celltype) %>% spread(study_celltype, logFC, fill=0)

correlation_matrix_escape <- cor(logFC_matrix_escape[,2:ncol(logFC_matrix_escape)], method = "spearman")

pdf('Aim_1/correlation_matrix_escape.pdf')
celltype = replace.names(gsub('_', '.', gsub('.+:', '', colnames(correlation_matrix_escape))))
column_ha = HeatmapAnnotation(disease = gsub(':.+', '', rownames(correlation_matrix_escape)),
celltype = celltype, col = list(disease=cols_study, celltype=cols_celltype))
row_ha = rowAnnotation(disease = gsub(':.+', '', colnames(correlation_matrix_escape)),
celltype = celltype, show_legend = c(FALSE, FALSE), col = list(disease=cols_study, celltype=cols_celltype))
Heatmap(correlation_matrix_escape, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), name='Rho',
show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha, right_annotation = row_ha)
dev.off()

### Jaccard index ###
jaccard_similarity <- function(set1, set2) {
  length(intersect(set1, set2)) / length(union(set1, set2))
}

MS_genes <- lapply(MS, function(x) x$gene)
names(MS_genes) <- paste('MS', names(MS_genes), sep=':')
pSS_genes <- lapply(pSS, function(x) x$gene)
names(pSS_genes) <- paste('pSS', names(pSS_genes), sep=':')
UC_genes <- lapply(UC, function(x) x$gene)
names(UC_genes) <- paste('UC', names(UC_genes), sep=':')
CD_colon_genes <- lapply(CD_colon, function(x) x$gene)
names(CD_colon_genes) <- paste('CD_colon', names(CD_colon_genes), sep=':')
CD_TI_genes <- lapply(CD_TI, function(x) x$gene)
names(CD_TI_genes) <- paste('CD_TI', names(CD_TI_genes), sep=':')
SLE_genes <- lapply(SLE, function(x) x$gene)
names(SLE_genes) <- paste('SLE', names(SLE_genes), sep=':')

# Combine all gene sets
all_genes <- c(MS_genes, pSS_genes, UC_genes, CD_colon_genes, CD_TI_genes, SLE_genes)
all_genes <- all_genes[gsub('.+:', '', names(all_genes)) %in% common_celltypes]

# Calculate Jaccard similarity for all combinations
results <- expand.grid(names(all_genes), names(all_genes))
colnames(results) <- c("CellType1", "CellType2")
results$JaccardIndex <- mapply(function(x, y) jaccard_similarity(all_genes[[x]], all_genes[[y]]),
                               results$CellType1, results$CellType2)
results <- spread(results, CellType2, JaccardIndex)
results[is.na(results)] <- 0

save(results, file='Aim_1/jaccard_similarity.RData')

pdf('Aim_1/jaccard_similarity.heatmap.pdf')
celltype = replace.names(gsub('_', '.', gsub('.+:', '', results$CellType1)))
column_ha = HeatmapAnnotation(disease = gsub(':.+', '', colnames(results)[-1]),
celltype = celltype, col = list(disease=cols_study, celltype=cols_celltype))
row_ha = rowAnnotation(disease = gsub(':.+', '', results$CellType1),
celltype = celltype, show_legend = c(FALSE, FALSE), col = list(disease=cols_study, celltype=cols_celltype))
Heatmap(as.matrix(results[,2:ncol(results)]), col=colorRamp2(c(0, 1), c("white", "red")), name='Jaccard Index',
show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha, right_annotation = row_ha)
dev.off()
