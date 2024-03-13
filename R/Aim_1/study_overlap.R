library(dplyr)
library(tidyr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)


source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')

MS <- deg.list('MS_GSE193770/differential.expression/edgeR', logfc=0.1)

pSS <- deg.list('pSS_GSE157278/differential.expression/edgeR', logfc=0.1)
UC <- deg.list('UC_GSE125527/differential.expression/edgeR', logfc=0.1)
CD_colon <- deg.list('CD_Kong/colon/differential.expression/edgeR', logfc=0.1)
names(CD_colon)[1] <- 'CD16-_NK_cells'
CD_TI <- deg.list('CD_Kong/TI/differential.expression/edgeR', logfc=0.1)
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
all_metrics$celltype[41] <- "CD16-.NK.cells"
all_metrics$celltype <- replace.names(gsub('_', '.', all_metrics$celltype))

write.table(all_metrics, 'Aim_1/study_metrics.txt', sep='\t', row.names=FALSE, quote=FALSE)

# Calculate enrichment of chrX in each cell type
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/fishers.test.degs.R')

MS_all <- deg.list('MS_GSE193770/differential.expression/edgeR', filter=FALSE)
pSS_all <- deg.list('pSS_GSE157278/differential.expression/edgeR', filter=FALSE)
UC_all <- deg.list('UC_GSE125527/differential.expression/edgeR', filter=FALSE)
CD_colon_all <- deg.list('CD_Kong/colon/differential.expression/edgeR', filter=FALSE)
CD_TI_all <- deg.list('CD_Kong/TI/differential.expression/edgeR', filter=FALSE)
SLE_all <- deg.list('lupus_Chun/differential.expression/edgeR', filter=FALSE)

fishers.all.chrX <- lapply(SLE_all, function(x) fisher.test.edgeR(x, rownames(chrX), 0.1, direction='none'))
fishers.all.chrX[sapply(fishers.all.chrX, function(x) x$p.value < 0.05)]

# Perform correlation of logFC between each combination of cell types
study <- pSS
correlation_matrix <- matrix(nrow = length(study), ncol = length(study),
                            dimnames = list(names(study), names(study)))

combinations <- combn(names(study), 2)

for(x in 1:ncol(combinations)) {
df1 <- study[[combinations[,x][1]]]
df2 <- study[[combinations[,x][2]]]

# Merge dataframes by gene
merged_df <- merge(df1, df2, by = "gene", suffixes = c(".1", ".2"))

# Calculate Spearman correlation
cor_result <- cor(merged_df$logFC.1, merged_df$logFC.2, method = "spearman")

# Save the result in the correlation matrix
correlation_matrix[combinations[,x][1], combinations[,x][2]] <- cor_result
correlation_matrix[combinations[,x][2], combinations[,x][1]] <- cor_result
}

diag(correlation_matrix) <- 1

# Change names
colnames(correlation_matrix) <- replace.names(gsub('_', '.', colnames(correlation_matrix)))
rownames(correlation_matrix) <- replace.names(gsub('_', '.', rownames(correlation_matrix)))

pdf(paste0('Aim_1/', 'pSS', '_correlation_matrix.pdf'))
Heatmap(correlation_matrix)
dev.off()

# Correlation of all cell types across all studies

MS_merged <- bind_rows(MS_all, .id = 'celltype')
MS_merged$study_celltype <- paste('MS', MS_merged$celltype, sep=':')
pSS_merged <- bind_rows(pSS_all, .id = 'celltype')
pSS_merged$study_celltype <- paste('pSS', pSS_merged$celltype, sep=':')
UC_merged <- bind_rows(UC_all, .id = 'celltype')
UC_merged$study_celltype <- paste('UC', UC_merged$celltype, sep=':')
CD_colon_merged <- bind_rows(CD_colon_all, .id = 'celltype')
CD_colon_merged$study_celltype <- paste('CD_colon', CD_colon_merged$celltype, sep=':')
CD_TI_merged <- bind_rows(CD_TI_all, .id = 'celltype')
CD_TI_merged$study_celltype <- paste('CD_TI', CD_TI_merged$celltype, sep=':')
SLE_merged <- bind_rows(SLE_all, .id = 'celltype')
SLE_merged$study_celltype <- paste('SLE', SLE_merged$celltype, sep=':')


all_merged <- bind_rows(MS_merged, pSS_merged, UC_merged, CD_colon_merged, CD_TI_merged, SLE_merged)

common_celltypes <- Reduce(intersect, list(names(pSS), names(UC), names(CD_colon), names(SLE)))
all_merged <- subset(all_merged, celltype %in% common_celltypes)

logFC_matrix <- all_merged[,c('gene', 'logFC', 'study_celltype')] %>% spread(study_celltype, logFC, fill=0)

correlation_matrix_all <- cor(logFC_matrix[,2:ncol(logFC_matrix)], method = "spearman")

library(RColorBrewer)
cols_study = brewer.pal(6, "Spectral")
names(cols_study) = c('CD_colon', 'CD_TI', 'MS', 'pSS', 'SLE', 'UC')
cols_celltype = brewer.pal(6, "Set3")
names(cols_celltype) = 

library(viridis)
cols_celltype = turbo(41)
names(cols_celltype) = unique(all_merged$celltype)

pdf('Aim_1/unfiltered_correlation_matrix.pdf')
column_ha = HeatmapAnnotation(disease = gsub(':.+', '', rownames(correlation_matrix_all)), 
celltype = gsub('.+:', '', colnames(correlation_matrix_all)), col = list(disease=cols_study, celltype=cols_celltype))
row_ha = rowAnnotation(disease = gsub(':.+', '', colnames(correlation_matrix_all)),
celltype = gsub('.+:', '', colnames(correlation_matrix_all)), show_legend = c(FALSE, FALSE), col = list(disease=cols_study, celltype=cols_celltype))
Heatmap(correlation_matrix_all, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), name='Rho',
show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha, right_annotation = row_ha)
dev.off()

logFC_matrix_sig <- all_merged %>% filter(FDR < 0.05 & abs(logFC) > 0.1) %>% select(gene, logFC, study_celltype) %>% spread(study_celltype, logFC, fill=0)
correlation_matrix_sig <- cor(logFC_matrix_sig[,2:ncol(logFC_matrix_sig)], method = "spearman")

pdf('Aim_1/sig_correlation_matrix.pdf')
column_ha = HeatmapAnnotation(disease = gsub(':.+', '', rownames(correlation_matrix_sig)), 
celltype = gsub('.+:', '', colnames(correlation_matrix_sig)), col = list(disease=cols_study, celltype=cols_celltype))
row_ha = rowAnnotation(disease = gsub(':.+', '', colnames(correlation_matrix_sig)),
celltype = gsub('.+:', '', colnames(correlation_matrix_sig)), show_legend = c(FALSE, FALSE), col = list(disease=cols_study, celltype=cols_celltype))
Heatmap(correlation_matrix_sig, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), name='Rho',
show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha, right_annotation = row_ha)
dev.off()


logFC_matrix_sig_chrX <- all_merged %>% filter(FDR < 0.05 & abs(logFC) > 0.1, gene %in% rownames(chrX)) %>% select(gene, logFC, study_celltype) %>% spread(study_celltype, logFC, fill=0)
correlation_matrix_sig_chrX <- cor(logFC_matrix_sig_chrX[,2:ncol(logFC_matrix_sig_chrX)], method = "spearman")

pdf('Aim_1/sig_correlation_matrix_chrX.pdf')
column_ha = HeatmapAnnotation(disease = gsub(':.+', '', rownames(correlation_matrix_sig_chrX)),
celltype = gsub('.+:', '', colnames(correlation_matrix_sig_chrX)), col = list(disease=cols_study, celltype=cols_celltype))
row_ha = rowAnnotation(disease = gsub(':.+', '', colnames(correlation_matrix_sig_chrX)),
celltype = gsub('.+:', '', colnames(correlation_matrix_sig_chrX)), show_legend = c(FALSE, FALSE), col = list(disease=cols_study, celltype=cols_celltype))
Heatmap(correlation_matrix_sig_chrX, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), name='Rho',
show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha, right_annotation = row_ha)
dev.off()

pdf('Aim_1/sig_logFC_matrix_chrX.pdf')
Heatmap(as.matrix(logFC_matrix_sig_chrX[,2:ncol(logFC_matrix_sig_chrX)]), name='logFC')
dev.off()

