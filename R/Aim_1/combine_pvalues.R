library(dplyr)
library(purrr)
library(metap)
library(metafor)
library(ComplexHeatmap)
library(circlize)

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')

chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')
tukiainen <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/Tukiainen.escape.txt')

pSS <- deg.list('pSS_GSE157278/differential.expression/edgeR/', filter=FALSE)
UC <- deg.list('UC_GSE125527/differential.expression/edgeR/', filter=FALSE)
CD_colon <- deg.list('CD_Kong/colon/differential.expression/edgeR/', filter=FALSE)
names(CD_colon)[2] <- 'CD16-_NK_cells'
CD_TI <- deg.list('CD_Kong/TI/differential.expression/edgeR/', filter=FALSE)
names(CD_TI)[1] <- 'CD16-_NK_cells'
SLE <- deg.list('lupus_Chun/differential.expression/edgeR/', filter=FALSE)
MS <- deg.list('MS_GSE193770/differential.expression/edgeR', filter=FALSE)

common <- Reduce(intersect, list(names(pSS), names(UC), names(CD_colon), names(CD_TI), names(SLE)))

combined_fdr_list <- list()
for(celltype in common){

    df_list <- list('pSS'=pSS[[celltype]], 'UC'=UC[[celltype]], 
        'CO'=CD_colon[[celltype]], 'TI'=CD_TI[[celltype]], 'SLE'=SLE[[celltype]],
        'MS'=MS[[celltype]])
    
    df_list <- df_list[sapply(df_list, function(x) is.null(x) == FALSE)]

    df <- df_list %>%
        imap(function(x, y) x %>% rename_with(~paste(., y, sep = '_'), -gene)) %>%
        reduce(full_join, by = 'gene') %>%
        data.frame()

    df_fdr <- df[,grep('FDR', colnames(df))]
    df_fdr[is.na(df_fdr)] <- 1

    df_logFC <- df[,grep('logFC', colnames(df))]
    df_logFC[is.na(df_logFC)] <- 0

    combined_fdr <- apply(df_fdr, 1, function(x) sumlog(x)$p)
    combined_logFC <- apply(df_logFC, 1, mean)

    combined_fdr_list[[celltype]] <- data.frame(gene = df$gene, 
        combined_fdr = combined_fdr, combined_logFC = combined_logFC)
}

save(combined_fdr_list, file = 'Aim_1/combined_fdr_list.Rdata')


degs_list <- lapply(combined_fdr_list, function(x){
    subset(x, combined_fdr < 0.01 & abs(combined_logFC) > 0.1)
})
degs <- bind_rows(degs_list, .id='celltype')

degs_wide <- dcast(gene ~ celltype, value=combined_logFC, data=degs)
degs_wide[is.na(degs_wide)] <- 0
rownames(degs_wide) <- degs_wide$gene
colnames(degs_wide) <- replace.names(gsub('_', '.', colnames(degs_wide)))
degs_wide <- degs_wide[,-1]


pdf('Aim_1/combined_logFC_sig_heatmap.pdf')
col <- colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red'))
Heatmap(as.matrix(degs_wide), col = col, name = 'combined_logFC',
show_row_names=FALSE) 
dev.off()



X.genes <- lapply(combined_fdr_list, function(x){
    x[x$gene %in% chrX & x$combined_fdr < 0.05,'gene']
})
load('../datasets/XCI/escapees.Rdata')
lapply(X.genes, function(x) length(x[x %in% rownames(escape)])/length(x))

lapply(X.genes, function(x) x[x %in% rownames(escape)])

fishers.test.combined <- function(data, genes, logfc){
    a <- nrow(data[data$combined_fdr < 0.05 & data$gene %in% genes & abs(data$combined_logFC) > logfc,])
    b <- nrow(data[data$combined_fdr < 0.05 & !data$gene %in% genes & abs(data$combined_logFC) > logfc,])
    c <- nrow(data[data$combined_fdr >= 0.05 & data$gene %in% genes,])
    d <- nrow(data[data$combined_fdr >= 0.05 & !data$gene %in% genes,])
    return(fisher.test(matrix(c(a, b, c, d), nrow = 2 ))$p.value)
}

lapply(combined_fdr_list, function(x) fishers.test.combined(x, rownames(escape), 0.1))

combined_fdr <- bind_rows(combined_fdr_list, .id = 'celltype')
subset_combined_fdr <- subset(combined_fdr, gene %in% chrX & combined_fdr < 0.05)
subset_combined_fdr$combined_logFC <- ifelse(subset_combined_fdr$combined_logFC == 'up', 1, -1)

plot_data <- reshape2::dcast(subset_combined_fdr, gene ~ celltype, value.var = 'combined_logFC')
plot_data[is.na(plot_data)] <- 0

library(ComplexHeatmap)
library(circlize)

pdf('meta_analysis_fdr.pdf')
col <- colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red'))
Heatmap(as.matrix(plot_data[,-1]), col = col, name = 'combined_logFC', row_names_gp = gpar(fontsize = 8))
dev.off()
    



df$combined_fdr <- apply(df[,-1], 1, function(x) sumlog(x)$p)


fishers.combined_fdr <- function(df, genes){
    a <- nrow(df[df$combined_fdr < 0.05 & df$gene %in% genes,])
    b <- nrow(df[df$combined_fdr < 0.05 & !df$gene %in% genes,])
    c <- nrow(df[df$combined_fdr >= 0.05 & df$gene %in% genes,])
    d <- nrow(df[df$combined_fdr >= 0.05 & !df$gene %in% genes,])
    return(fisher.test(matrix(c(a, b, c, d), nrow = 2 ))$p.value)
}



df_subset <- subset(df, combined_fdr < 0.05)

df_subset[df_subset$gene %in% chrX,]