library(dplyr)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(reshape2)
library(UpSetR)
library(org.Hs.eg.db)
library(fgsea)

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/wFisher.R')

chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')
tukiainen <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/Tukiainen.escape.txt')
X.immune <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/X.immune.txt')

# Read in propellor results from each study
pSS_prop <- readRDS('pSS_GSE157278/propellor.asin.RDS')
UC_prop <- readRDS('UC_GSE125527/propellor.asin.RDS')
CO_prop <- readRDS('CD_Kong/colon/propeller.asin.RDS')
TI_prop <- readRDS('CD_Kong/TI/propeller.asin.RDS')
SLE_prop <- read.delim('lupus_Chun/propeller.asin.condition.abundance.txt')
rownames(SLE_prop) <- SLE_prop[,1]
MS_prop <- readRDS('MS_GSE193770/propellor.asin.RDS')

prop_list <- list(pSS_prop, UC_prop, CO_prop, TI_prop, SLE_prop, MS_prop)

save(prop_list, file = 'Aim_1/prop_list.Rdata')

all_celltypes <- unique(unlist(lapply(prop_list, rownames)))
propellor_matrix <- matrix(0, nrow = length(all_celltypes), ncol = 6)
rownames(propellor_matrix) <- all_celltypes
colnames(propellor_matrix) <- c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS')
# Add Tstatistic to matrix
for(i in 1:6){
    propellor_matrix[which(rownames(propellor_matrix) %in% rownames(prop_list[[i]])), i] <- prop_list[[i]][,'Tstatistic']
}

propellor_fdr <- matrix(1, nrow = length(all_celltypes), ncol = 6)
rownames(propellor_fdr) <- all_celltypes
colnames(propellor_fdr) <- c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS')
# Add FDR to matrix
for(i in 1:6){
    propellor_fdr[which(rownames(propellor_fdr) %in% rownames(prop_list[[i]])), i] <- prop_list[[i]][,'FDR']
}
significance_matrix <- ifelse(propellor_fdr < 0.05, "*", "ns")
significance_matrix <- ifelse(propellor_fdr < 0.01, "**", significance_matrix)
significance_matrix <- ifelse(propellor_fdr < 0.001, "***", significance_matrix)

# Heatmap of T-statistic with FDR astrix overlay
pdf('Aim_1/propellor_heatmap.pdf')
col <- colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red'))
Heatmap(propellor_matrix, col = col, name = 'T-statistic', show_row_names=TRUE, show_column_names=TRUE,
row_names_gp = gpar(fontsize = 8),
cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(significance_matrix[i, j], x, y, gp = gpar(fontsize = 5))
        })
dev.off()


# Read in edgeR results for each study
pSS <- deg.list('pSS_GSE157278/differential.expression/edgeR/', filter=FALSE)
UC <- deg.list('UC_GSE125527/differential.expression/edgeR/', filter=FALSE)
CO <- deg.list('CD_Kong/colon/differential.expression/edgeR/', filter=FALSE)
names(CO)[2] <- 'CD16-_NK_cells'
TI <- deg.list('CD_Kong/TI/differential.expression/edgeR/', filter=FALSE)
names(TI)[1] <- 'CD16-_NK_cells'
SLE <- deg.list('lupus_Chun/differential.expression/edgeR/', filter=FALSE)
MS <- deg.list('MS_GSE193770/differential.expression/edgeR', filter=FALSE)

# Add studies to list
study_list <- list(pSS, UC, CO, TI, SLE, MS)

study_celltypes <- lapply(study_list, names)
names(study_celltypes) <- c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS')

celltypes_mtx <- fromList(study_celltypes)
rownames(celltypes_mtx) <- unique(unlist(study_celltypes))
### select common cell types with at least 5 samples ###
common <- rownames(celltypes_mtx)[rowSums(celltypes_mtx) >= 5]

# Combine p-values and logFC with wFisher
sample_size <- list('pSS'=10, 'UC'=8, 'CO'=13, 'TI'=16, 'SLE'=229, 'MS'=6)

combined_fdr_list <- list()
for(celltype in common){

    df_list <- list('pSS'=pSS[[celltype]], 'UC'=UC[[celltype]], 
        'CO'=CO[[celltype]], 'TI'=TI[[celltype]], 'SLE'=SLE[[celltype]],
        'MS'=MS[[celltype]])
    
    df_list <- df_list[sapply(df_list, function(x) is.null(x) == FALSE)]

    df <- df_list %>%
        imap(function(x, y) x %>% rename_with(~paste(., y, sep = '_'), -gene)) %>%
        purrr::reduce(full_join, by = 'gene') %>%
        data.frame()

    df[,grep('FDR', colnames(df))][is.na(df[,grep('FDR', colnames(df))])] <- 1

    df[,grep('logFC', colnames(df))][is.na(df[,grep('logFC', colnames(df))])] <- 0
    df[,grep('logFC', colnames(df))] <- ifelse(df[,grep('logFC', colnames(df))] > 0.1, 1, 
                                        ifelse(df[,grep('logFC', colnames(df))] < -0.1, -1, 0))

    result <- apply(df, 1, function(x){
        wFisher(p=as.numeric(x[grep('FDR', names(x))]), 
        weight = unlist(sample_size[names(df_list)]), 
        is.onetail=FALSE, 
        eff.sign = x[grep('logFC', names(x))])
    })
    combined_fdr <- result[[1]]
    combined_logFC <- apply(df[,grep('logFC', colnames(df))], 1, function(x) mean(x))

    combined_fdr_list[[celltype]] <- data.frame(gene = df$gene, 
        combined_fdr = unlist(lapply(result, function(x) x$p)),
        combined_logFC = unlist(lapply(result, function(x) x$overall.eff.direction)))
}

save(combined_fdr_list, file = 'Aim_1/combined_fdr_list.Rdata')
load('Aim_1/combined_fdr_list.Rdata')

# Add lists together
fdr_df <- bind_rows(combined_fdr_list, .id='celltype')
# Reshape data
fdr_df_wide <- dcast(gene ~ celltype, value.var='combined_fdr', data=fdr_df)
rownames(fdr_df_wide) <- fdr_df_wide$gene
fdr_df_wide <- fdr_df_wide[,-1]
colnames(fdr_df_wide) <- replace.names(gsub('_', '.', colnames(fdr_df_wide)))
# Replace NA with 1
fdr_df_wide[is.na(fdr_df_wide)] <- 1

# Remove genes not significant in any cell type
fdr_df_wide <- fdr_df_wide[apply(fdr_df_wide, 1, function(x) length(which(x < 0.05)) > 1),]

pdf('Aim_1/combined_fdr_heatmap.pdf')
col <- colorRamp2(c(0, 0.05, 1), c('red', 'blue', 'white'))
Heatmap(as.matrix(fdr_df_wide), col = col, name = 'FDR', show_row_names=FALSE, use_raster=TRUE)
dev.off()


# Test for enrichment of escape genes
chisq.test.combined <- function(data, genes){
    a <- nrow(data[data$combined_fdr < 0.05 & data$gene %in% genes,])
    b <- nrow(data[data$combined_fdr < 0.05 & !data$gene %in% genes,])
    c <- nrow(data[data$combined_fdr >= 0.05 & data$gene %in% genes,])
    d <- nrow(data[data$combined_fdr >= 0.05 & !data$gene %in% genes,])
    return(chisq.test(matrix(c(a, b, c, d), nrow = 2 ))$p.value)
}

### XCI escape ###
deg_escape <- lapply(combined_fdr_list, function(x){
    data.frame(
        enrichment = chisq.test.combined(x, rownames(escape)),
        size = nrow(x[x$combined_fdr < 0.05 & x$gene %in% rownames(escape),])
    )
})
deg_escape <- bind_rows(deg_escape, .id='celltype')
deg_escape$celltype <- replace.names(gsub('_', '.', deg_escape$celltype))

pdf('Aim_1/escape_enrichment.pdf')
ggplot(deg_escape, aes(x=celltype, y=-log10(enrichment), size=size, color=-log10(enrichment))) +
    geom_point() +
    scale_color_gradient(low='blue', high='red') +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    theme_minimal() +
    xlab('') +
    ggtitle('XCI escape genes') +
    coord_flip()
dev.off()

### Interferome ### 
dir.create('Aim_1/combined_degs')
lapply(names(combined_fdr_list), function(x){
    tmp <- subset(combined_fdr_list[[x]], combined_fdr < 0.05)
    write.csv(tmp, file = paste0('Aim_1/combined_degs/', x, '.csv'), row.names=FALSE)
})

interferome <- list.files('Aim_1/interferome', full.names=TRUE) %>%
    lapply(function(x) read.delim(x, skip=19))
names(interferome) <- gsub('.*/|\\.txt', '', list.files('Aim_1/interferome'))

deg_interferome <- lapply(names(combined_fdr_list), function(x){
    data.frame(
        enrichment = chisq.test.combined(combined_fdr_list[[x]], interferome[[x]]$Gene.Name),
        size = nrow(subset(combined_fdr_list[[x]], combined_fdr < 0.05 & gene %in% interferome[[x]]$Gene.Name))
    )
})
names(deg_interferome) <- replace.names(gsub('_', '.', names(combined_fdr_list)))
deg_interferome <- bind_rows(deg_interferome, .id='celltype')

pdf('Aim_1/interferome_enrichment.pdf')
ggplot(deg_interferome, aes(x=celltype, y=-log10(enrichment), size=size, color=-log10(enrichment))) +
    geom_point() +
    scale_color_gradient(low='blue', high='red') +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    theme_minimal() +
    xlab('') +
    ggtitle('Interferome') +
    coord_flip()
dev.off()

### Heatmap of XCI escape degs ###
XCI_deg <- lapply(combined_fdr_list, function(x){
    subset(x, combined_fdr < 0.05 & gene %in% rownames(escape))
})
XCI_deg <- bind_rows(XCI_deg, .id='celltype')
XCI_deg$combined_logFC <- ifelse(XCI_deg$combined_logFC == '-', 1, 0)
XCI_deg$celltype <- replace.names(gsub('_', '.', XCI_deg$celltype))
XCI_deg_wide <- dcast(gene ~ celltype, value.var='combined_logFC', data=XCI_deg)
XCI_deg_wide[is.na(XCI_deg_wide)] <- 0

pdf('Aim_1/XCI_deg_heatmap.pdf')
rownames(XCI_deg_wide) <- XCI_deg_wide$gene
ha = rowAnnotation(escape = anno_mark(at = which(rownames(XCI_deg_wide) %in% X.immune[,1]), 
    labels = rownames(XCI_deg_wide)[rownames(XCI_deg_wide) %in% X.immune[,1]]))
col <- colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red'))
Heatmap(as.matrix(XCI_deg_wide[,-1]), col = col, name = 'logFC', show_row_names=FALSE, right_annotation = ha)
dev.off()


# Map gene to chromosome
combined_fdr_list_sig <- lapply(combined_fdr_list, function(x){
    tmp <- subset(x, combined_fdr < 0.05)
    tmp$chrom <- mapIds(org.Hs.eg.db, keys = tmp$gene, column = "CHR", keytype = "SYMBOL", multiVals = "first")
    tmp <- subset(tmp, chrom != NA)
    tmp$chrom <- factor(tmp$chrom, levels=c(1:22, 'X', 'Y'))
    return(tmp)
})
head(combined_fdr_list_sig[[1]])
combined_fdr_sig <- bind_rows(combined_fdr_list_sig, .id='celltype')

fdr_df_sig <- subset(fdr_df, combined_fdr < 0.05)
chrom_info <- mapIds(org.Hs.eg.db, keys = unique(fdr_df$gene), column = "CHR", keytype = "SYMBOL", multiVals = "first")
fdr_df_sig$chrom <- as.character(chrom_info[fdr_df_sig$gene])
fdr_df_sig <- fdr_df_sig[which(fdr_df_sig$chrom != 'NA'),]
fdr_df_sig$chrom <- factor(fdr_df_sig$chrom, levels=c(1:22, 'X', 'Y'))
fdr_df_sig$celltype <- replace.names(gsub('_', '.', fdr_df_sig$celltype))

# Remove Y chromosome
fdr_df_sig <- subset(fdr_df_sig, chrom != 'Y')

save(fdr_df_sig, file = 'Aim_1/fdr_df_sig.Rdata')

pdf('Aim_1/chromosome_freq_barplot.pdf')
ggplot(fdr_df_sig, aes(x=chrom)) +
    geom_bar() +
    theme_minimal() +
    xlab('Chromosome') +
    ylab('Frequency') + 
    theme(axis.text.x = element_text(size=5)) +
    facet_wrap(~celltype, ncol=3, nrow=4, scales='free_y', strip.position = "bottom") 
dev.off()

### Gene set enrichment analysis
### GSEA hallmark ###
hallmark <- gmtPathways('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/h.all.v7.5.1.symbols.gmt')
position <- gmtPathways('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/c1.all.v2023.1.Hs.symbols.gmt')

load('pSS_GSE157278/Aim_1_2024/figure.data/fgsea_list.RData')
pSS_fgsea <- bind_rows(fgsea_list, .id='celltype')
load('UC_GSE125527//Aim_1_2024/figure.data/fgsea_list.RData')
UC_fgsea <- bind_rows(fgsea_list, .id='celltype')
load('CD_Kong/colon/Aim_1_2024/figure.data/fgsea_list.RData')
CO_fgsea <- bind_rows(fgsea_list, .id='celltype')
load('CD_Kong/TI/Aim_1_2024/figure.data/fgsea_list.RData')
TI_fgsea <- bind_rows(fgsea_list, .id='celltype')
load('lupus_Chun/Aim_1_2024/figure.data/fgsea_list.RData')
SLE_fgsea <- bind_rows(fgsea_list, .id='celltype')
load('MS_GSE193770/Aim_1_2024/figure.data/fgsea_list.RData')
MS_fgsea <- bind_rows(fgsea_list, .id='celltype')

fgsea_list <- list('pSS'=pSS_fgsea, 'UC'=UC_fgsea, 'CO'=CO_fgsea, 'TI'=TI_fgsea, 'SLE'=SLE_fgsea, 'MS'=MS_fgsea)
save(fgsea_list, file='Aim_1/fgsea_list.RData')

fgsea_df <- bind_rows(fgsea_list, .id='study')

fgsea_df <- subset(fgsea_df, celltype %in% replace.names(gsub('_', '.', names(combined_fdr_list))))
fgsea_df$celltype_study <- paste(fgsea_df$celltype, fgsea_df$study, sep=':')

fgsea_wide <- dcast(pathway ~ celltype_study, value.var='NES', data=fgsea_df)
rownames(fgsea_wide) <- gsub('HALLMARK_', '', fgsea_wide$pathway)
fgsea_wide <- fgsea_wide[,-1]
fgsea_wide[is.na(fgsea_wide)] <- 0
study <- gsub('.+:', '', colnames(fgsea_wide))
colnames(fgsea_wide) <- gsub(':.+', '', colnames(fgsea_wide))

pdf('Aim_1/fgsea_heatmap.pdf')
anno <- HeatmapAnnotation(study = study, col = list(study = c('pSS'='#E69F00',  # orange
                                                              'UC'='#009E73',  # light blue
                                                              'CO'='#56B4E9',  # teal
                                                              'TI'='#0072B2',  # darker teal
                                                              'SLE'='#D55E00',  # reddish orange
                                                              'MS'='#CC79A7')))
col <- colorRamp2(c(-2, 0, 2), c('#67a9cf', 'white', '#ef8a62'))  # light blue to white to coral red
Heatmap(as.matrix(fgsea_wide), col = col, name = 'NES', show_row_names=TRUE, show_column_names=TRUE,
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 5),
        top_annotation = anno)
dev.off()
