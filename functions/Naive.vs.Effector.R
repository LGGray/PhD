library(edgeR)
library(Seurat)
library(MAST)
library(qvalue)
library(tidyverse)

if(dir.exists('differential.expression/naive_effector') != TRUE){dir.create('differential.expression/naive_effector')}

cell_1 <- commandArgs(trailingOnly=TRUE)[1]
cell_2 <- commandArgs(trailingOnly=TRUE)[2]

pbmc <- readRDS("pbmc.female.RDS")

# Control comparison
pbmc.control <- subset(pbmc, cellTypist %in% c(cell_1, cell_2) & condition == 'control')
keep <- rowSums(pbmc.control@assays$RNA@counts > 0) > ncol(pbmc.control) * 0.05
features <- names(keep[keep == T])
pbmc.control <- subset(pbmc.control, features=features)
# MAST differential expression
control <- FindMarkers(pbmc.control, slot='counts', ident.1 = cell_1, ident.2 = cell_2,
                             test.use = "MAST", latent.vars='SV1', min.pct = 0, logfc.threshold = 0)
control <- cbind(gene = rownames(control), control)

# Disease comparison
pbmc.disease <- subset(pbmc, cellTypist %in% c(cell_1, cell_2) & condition == 'disease')
keep <- rowSums(pbmc.disease@assays$RNA@counts > 0) > ncol(pbmc.disease) * 0.05
features <- names(keep[keep == T])
pbmc.disease <- subset(pbmc.disease, features=features)
# MAST differential expression
disease <- FindMarkers(pbmc.disease, slot='counts', ident.1 = cell_1, ident.2 = cell_2,
                             test.use = "MAST", latent.vars='SV1', min.pct = 0, logfc.threshold = 0)
disease <- cbind(gene = rownames(disease), disease)


merged <- merge(control, disease, by='gene', suffixes = c('.control', '.disease'))

cor(merged$avg_log2FC.control, merged$avg_log2FC.disease)
cor(merged$p_val_adj.control, merged$p_val_adj.disease)

merged.deg <- subset(merged, p_val_adj.control < 0.05 & p_val_adj.disease < 0.05)
cor(merged.deg$avg_log2FC.control, merged.deg$avg_log2FC.disease)
cor(merged.deg$p_val_adj.control, merged.deg$p_val_adj.disease)

control.chrX <- subset(control, p_val_adj < 0.05 & gene %in% rownames(chrX))
nrow(control.chrX)
disease.chrX <- subset(disease, p_val_adj < 0.05 & gene %in% rownames(chrX))
nrow(disease.chrX)

chisq.test.MAST(control, rownames(chrX), 0)
chisq.test.MAST(disease, rownames(chrX), 0)