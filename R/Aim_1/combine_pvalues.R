library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(UpSetR)
library(org.Hs.eg.db)
library(fgsea)
library(scico)

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
pSS_prop <- readRDS('pSS_GSE157278/propeller.asin.RDS')
UC_prop <- readRDS('UC_GSE125527/propeller.asin.RDS')
CO_prop <- readRDS('CD_Kong/colon/propeller.asin.RDS')
TI_prop <- readRDS('CD_Kong/TI/propeller.asin.RDS')
SLE_prop <- readRDS('lupus_Chun/propeller.asin.RDS')
MS_prop <- readRDS('MS_GSE193770/propeller.asin.RDS')

prop_list <- list(pSS_prop, UC_prop, CO_prop, TI_prop, SLE_prop, MS_prop)

save(prop_list, file = 'Aim_1/prop_list.Rdata')

extract_tstat <- function(df, name) {
  df <- df[, c("BaselineProp.clusters", "Tstatistic")]
  names(df)[2] <- name
  return(df)
}
tstat_list <- Map(extract_tstat, prop_list, 1:6)
propeller_matrix <- Reduce(function(x, y) merge(x, y, by = "BaselineProp.clusters", all = TRUE), tstat_list)
rownames(propeller_matrix) <- propeller_matrix$BaselineProp.clusters
propeller_matrix <- propeller_matrix[,-1]
colnames(propeller_matrix) <- c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS')
propeller_matrix[is.na(propeller_matrix)] <- 0
propeller_matrix <- as.matrix(propeller_matrix)

extract_FDR <- function(df, name) {
  df <- df[, c("BaselineProp.clusters", "FDR")]
  names(df)[2] <- name
  return(df)
}
fdr_list <- Map(extract_FDR, prop_list, 1:6)
propeller_fdr <- Reduce(function(x, y) merge(x, y, by = "BaselineProp.clusters", all = TRUE), fdr_list)
rownames(propeller_fdr) <- propeller_fdr$BaselineProp.clusters
propeller_fdr <- propeller_fdr[,-1]
colnames(propeller_fdr) <- c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS')
propeller_fdr <- ifelse(propeller_fdr < 0.05, "*", "ns")
propeller_fdr <- ifelse(propeller_fdr < 0.01, "**", propeller_fdr)
propeller_fdr <- ifelse(propeller_fdr < 0.001, "***", propeller_fdr)

# Heatmap of T-statistic with FDR astrix overlay
pdf('Aim_1/propellor_heatmap.pdf')
col <- colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red'))
Heatmap(propeller_matrix, col = col, name = 'T-statistic', show_row_names=TRUE, show_column_names=TRUE,
row_names_gp = gpar(fontsize = 8),
cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(propeller_fdr[i, j], x, y, gp = gpar(fontsize = 5))
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
study_list <- list('pSS'=pSS, 'UC'=UC, 'CO'=CO, 'TI'=TI, 'SLE'=SLE, 'MS'=MS)

study_celltypes <- lapply(study_list, names)
names(study_celltypes) <- c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS')
celltypes_mtx <- fromList(study_celltypes)
rownames(celltypes_mtx) <- unique(unlist(study_celltypes))
### select common cell types with at least 5 samples ###
common <- rownames(celltypes_mtx)[rowSums(celltypes_mtx) >= 5]

### Identify Jaccard index of upregulated genes between common cell types ###
jaccard_index <- function(set1, set2){
    set1 <- set1[!is.na(set1)]
    set2 <- set2[!is.na(set2)]
    intersect <- length(intersect(set1, set2))
    union <- length(union(set1, set2))
    return(intersect/union)
}

jaccard_matrix_list <- list()
for(i in 1:length(common)){
    celltype <- lapply(study_list, function(x) x[[common[i]]])
    degs <- lapply(celltype, function(x){
        if(is.null(x)){
            return(NULL)
        }
        subset(x, abs(logFC) > 0.1 & FDR < 0.05)$gene
    })

    jaccard_matrix <- matrix(NA, nrow=length(study_list), ncol=length(study_list))
    rownames(jaccard_matrix) <- names(study_list)
    colnames(jaccard_matrix) <- names(study_list)

    for(n in 1:length(study_list)){
        for(m in 1:length(study_list)){
            jaccard_matrix[n,m] <- jaccard_index(degs[[n]], degs[[m]])
        }
    }
    jaccard_matrix_list[[common[i]]] <- jaccard_matrix
}

heatmaps <- lapply(names(jaccard_matrix_list), function(x){
    Heatmap(jaccard_matrix_list[[x]], name='Jaccard index', show_row_names=TRUE, show_column_names=TRUE,
    cluster_columns=FALSE, cluster_rows=FALSE, column_title=x, column_title_gp = gpar(fontsize = 5))
})
pdf('Aim_1/jaccard_heatmaps.pdf', width=20, height=5)
heatmaps[[1]] + heatmaps[[2]] + heatmaps[[3]] + heatmaps[[4]] + heatmaps[[5]] + heatmaps[[6]] + heatmaps[[7]] + heatmaps[[8]] + heatmaps[[9]] + heatmaps[[10]] + heatmaps[[11]]
dev.off()


### Combine p-values and logFC with wFisher
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
# Reshape the data frame
fdr_df_wide <- dcast(gene ~ celltype, value.var='combined_fdr', data=fdr_df)
rownames(fdr_df_wide) <- fdr_df_wide$gene
fdr_df_wide <- as.matrix(fdr_df_wide[,-1])
colnames(fdr_df_wide) <- replace.names(gsub('_', '.', colnames(fdr_df_wide)))
# Replace NA with 1
fdr_df_wide[is.na(fdr_df_wide)] <- 1

# Remove genes not significant in any cell type
fdr_df_wide <- fdr_df_wide[apply(fdr_df_wide, 1, function(x) length(which(x < 0.05)) > 1),]

pdf('Aim_1/combined_fdr_heatmap.pdf')
col <- colorRamp2(c(0, 0.05, 1), c('red', 'blue', 'white'))
Heatmap(as.matrix(fdr_df_wide), col = col, name = 'FDR', show_row_names=FALSE, use_raster=TRUE)
dev.off()

# Upset plot of significant genes across cell types
degs <- lapply(combined_fdr_list, function(x){
    subset(x, combined_fdr < 0.05)$gene
})
names(degs) <- replace.names(gsub('_', '.', names(combined_fdr_list)))
degs_mtx <- fromList(degs)
rownames(degs_mtx) <- unique(unlist(degs))
pdf('Aim_1/upset_plot.pdf', onefile=F, width=10, height=10)
upset(degs_mtx, order.by = 'freq', nsets = 11)
dev.off()

# Test for enrichment of gene set
chisq.test.combined <- function(data, genes){
    a <- nrow(data[data$combined_fdr < 0.05 & data$gene %in% genes,])
    b <- nrow(data[data$combined_fdr < 0.05 & !data$gene %in% genes,])
    c <- nrow(data[data$combined_fdr >= 0.05 & data$gene %in% genes,])
    d <- nrow(data[data$combined_fdr >= 0.05 & !data$gene %in% genes,])
    return(chisq.test(matrix(c(a, b, c, d), nrow = 2 ))$p.value)
}

### X chromosome genes ###
deg_X <- lapply(combined_fdr_list, function(x){
    data.frame(
        enrichment = chisq.test.combined(x, chrX),
        size = nrow(x[x$combined_fdr < 0.05 & x$gene %in% chrX,])
    )
})
deg_X <- bind_rows(deg_X, .id='celltype')
deg_X$celltype <- replace.names(gsub('_', '.', deg_X$celltype))

pdf('Aim_1/chrX_enrichment.pdf')
ggplot(deg_X, aes(x=celltype, y=-log10(enrichment), size=size, color=-log10(enrichment))) +
    geom_point() +
    scale_color_gradient(low='blue', high='red') +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    theme_minimal() +
    xlab('') +
    ggtitle('Enrichment of X chromosome genes') +
    coord_flip()
dev.off()

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
    ggtitle('Enrichment of XCI escape genes') +
    coord_flip()
dev.off()

### Set colours and shape for celltypes and study
celltypes_colour <- scico(length(common), palette = 'lipari')
celltypes_colour <- setNames(celltypes_colour, replace.names(gsub('_', '.', common)))
study_colours <- scico(6, palette = 'roma')
study_colours <- setNames(study_colours, c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS'))
study_shapes <- c('pSS'=21, 'UC'=22, 'CO'=23, 'TI'=23, 'SLE'=24, 'MS'=25)

### Plot common features ###
features <- c('TSIX', 'FTX', 'GABRE', 'EIF2S3', 'SAT1', 'TSC22D3', 'OGT', 'KDM6A')
deg_escape <- lapply(combined_fdr_list, function(x){
    subset(x, combined_fdr < 0.05 & gene %in% rownames(escape))$gene
})
feature_celltypes <- names(deg_escape[sapply(deg_escape, function(x) features %in% x)])

features_df <- lapply(study_list, function(study){
    tmp <- lapply(study, function(celltype){
        subset(celltype, gene %in% features)[,c('gene', 'logFC', 'FDR')]
    })
    df <- bind_rows(tmp, .id='celltype')
    return(df)
})
names(features_df) <- c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS')
features_df <- bind_rows(features_df, .id='study')
plot.data <- subset(features_df, gene %in% features & celltype %in% feature_celltypes)
plot.data$celltype <- factor(replace.names(gsub('_', '.', plot.data$celltype)))
plot.data$study <- factor(plot.data$study)
plot.data$gene <- factor(plot.data$gene, levels=c('TSIX', 'FTX', 'GABRE', 'EIF2S3', 'SAT1', 'TSC22D3', 'OGT', 'KDM6A'))

pdf('Aim_1/common_features_volcano.pdf', width=10, height=10)
ggplot(plot.data, aes(x=logFC, y=-log10(FDR), color=celltype, shape=study)) +
    geom_point() +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    geom_vline(xintercept=c(-0.1, 0.1), linetype='dashed') +
    scale_color_manual(values=celltypes_colour) +
    scale_shape_manual(values=study_shapes) +
    facet_wrap(~gene, scales='free_x', ncol=4)
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

# Find cell type specific genes ##
deg_escape <- lapply(combined_fdr_list, function(x){
    subset(x, combined_fdr < 0.05 & gene %in% rownames(escape))$gene
})
deg_escape_mtx <- fromList(deg_escape)
rownames(deg_escape_mtx) <- unique(unlist(deg_escape))
rownames(deg_escape_mtx[rowSums(deg_escape_mtx) == 1,])

### XIST RBP ###
XIST_RBP <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/XIST_RBP.txt')
XIST_RBP$Gene.Symbol

library(dplyr)
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

convert_mouse_to_human <- function(gene_list){

  output = c()

  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }

  return (output)
}

XIST_RBP <- convert_mouse_to_human(XIST_RBP$Gene.Symbol)
deg_XIST_RBP <- lapply(combined_fdr_list, function(x){
    data.frame(
        enrichment = chisq.test.combined(x, XIST_RBP),
        size = nrow(x[x$combined_fdr < 0.05 & x$gene %in% XIST_RBP,])
    )
})
deg_XIST_RBP <- bind_rows(deg_XIST_RBP, .id='celltype')
deg_XIST_RBP$celltype <- replace.names(gsub('_', '.', deg_XIST_RBP$celltype))

pdf('Aim_1/XIST_RBP_enrichment.pdf')
ggplot(deg_XIST_RBP, aes(x=celltype, y=-log10(enrichment), size=size, color=-log10(enrichment))) +
    geom_point() +
    scale_color_gradient(low='blue', high='red') +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    theme_minimal() +
    xlab('') +
    ggtitle('XIST RNA Binding Proteins') +
    coord_flip()
dev.off()

deg_XIST_RBP <- lapply(combined_fdr_list, function(x){
    subset(x, combined_fdr < 0.05 & gene %in% XIST_RBP)$gene
})
sort(table(unlist(deg_XIST_RBP)))

deg_XIST_RBP[sapply(deg_XIST_RBP, function(x) 'PRC2' %in% x)]

# DisGeneNet
disgenet_pSS <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/pSS.tsv')$Gene
disgenet_UC <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/UC.tsv')$Gene
disgenet_CO <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/CD_colon.tsv')$Gene
disgenet_TI <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/CD_TI.tsv')$Gene
disgenet_SLE <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/SLE.tsv')$Gene
disgenet_MS <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/MS.tsv')$Gene
disgenet_list <- list('pSS'=disgenet_pSS, 'UC'=disgenet_UC, 'CO'=disgenet_CO, 'TI'=disgenet_TI, 'SLE'=disgenet_SLE, 'MS'=disgenet_MS)

disgenet <- lapply(disgenet_list, function(x){
    tmp <- lapply(combined_fdr_list, function(y){
        data.frame(
            enrichment=chisq.test.combined(y, x),
            size=nrow(y[y$combined_fdr < 0.05 & y$gene %in% x,]))
    })
    bind_rows(tmp, .id='celltype')
})
disgenet_df <- bind_rows(disgenet, .id='study')
disgenet_df$celltype <- replace.names(gsub('_', '.', disgenet_df$celltype))
disgenet$study <- factor(disgenet_df$study, levels=c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS'))

pdf('Aim_1/disgenet_enrichment.pdf')
ggplot(disgenet_df, aes(x=celltype, y=-log10(enrichment), size=size, color=-log10(enrichment))) +
    geom_point() +
    scale_color_gradient(low='blue', high='red') +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    theme_minimal() +
    xlab('') +
    ggtitle('DisGeNet') +
    coord_flip() +
    facet_wrap(~study, ncol=3)
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

fdr_df <- bind_rows(combined_fdr_list, .id='celltype')
fdr_df_sig <- subset(fdr_df, combined_fdr < 0.05)
chrom_info <- mapIds(org.Hs.eg.db, keys = unique(fdr_df$gene), column = "CHR", keytype = "SYMBOL", multiVals = "first")
fdr_df_sig$chrom <- as.character(chrom_info[fdr_df_sig$gene])
fdr_df_sig <- fdr_df_sig[which(fdr_df_sig$chrom != 'NA'),]
fdr_df_sig$chrom <- factor(fdr_df_sig$chrom, levels=c(1:22, 'X', 'Y'))
fdr_df_sig$celltype <- replace.names(gsub('_', '.', fdr_df_sig$celltype))

# Remove Y chromosome
fdr_df_sig <- subset(fdr_df_sig, chrom != 'Y')

range(unlist(lapply(split(fdr_df_sig, fdr_df_sig$celltype), function(x){
    which(names(sort(table(x$chrom), decreasing=TRUE))=='X')
})))

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
load('Aim_1/fgsea_list.RData')

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

foo <- subset(fgsea_df, pathway %in% 'HALLMARK_OXIDATIVE_PHOSPHORYLATION')
foo[order(foo$celltype_study),]
unique(foo$celltype)
unique(foo$study)