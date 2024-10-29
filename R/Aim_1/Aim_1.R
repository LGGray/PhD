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
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/chisq.test.degs.R')

chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')
X.immune <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/X.immune.txt')
X.gwas <- read.csv('/directflow/SCCGGroupShare/projects/lacgra/gwas/chrX.variants.csv')
X.gwas <- unique(unlist(strsplit(X.gwas$mappedGenes, ',')))

# Read in propellor results from each study
pSS_prop <- readRDS('pSS_GSE157278/propeller.asin.RDS')
UC_prop <- readRDS('UC_GSE182270/propeller.asin.RDS')
CO_prop <- readRDS('CD_Kong/colon/propeller.asin.RDS')
TI_prop <- readRDS('CD_Kong/TI/propeller.asin.RDS')
SLE_prop <- readRDS('lupus_Chun/propeller.asin.RDS')
MS_prop <- readRDS('MS_GSE193770/propeller.asin.RDS')

prop_list <- list(pSS_prop, UC_prop, CO_prop, TI_prop, SLE_prop, MS_prop)

lapply(prop_list, function(x) {
  paste(nrow(x[x$FDR < 0.05,]), nrow(x))
})

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
UC <- deg.list('UC_GSE182270/differential.expression/edgeR/', filter=FALSE)
names(UC)[3] <- 'CD16-_NK_cells'
CO <- deg.list('CD_Kong/colon/differential.expression/edgeR/', filter=FALSE)
names(CO)[2] <- 'CD16-_NK_cells'
TI <- deg.list('CD_Kong/TI/differential.expression/edgeR/', filter=FALSE)
names(TI)[1] <- 'CD16-_NK_cells'
SLE <- deg.list('lupus_Chun/differential.expression/edgeR/', filter=FALSE)
MS <- deg.list('MS_GSE193770/differential.expression/edgeR', filter=FALSE)

# Add studies to list
study_list <- list('pSS'=pSS, 'UC'=UC, 'CO'=CO, 'TI'=TI, 'SLE'=SLE, 'MS'=MS)

# Save the list
save(study_list, file = 'Aim_1/study_list.Rdata')

# All cell types
all_celltypes <- replace.names(gsub('_', '.', unique(unlist(lapply(study_list, function(x) names(x))))))

### Set colours and shape for celltypes and study
celltypes_colour <- scico(length(all_celltypes), palette = 'lipari')
celltypes_colour <- setNames(celltypes_colour, all_celltypes)
study_colours <- scico(6, palette = 'roma')
study_colours <- setNames(study_colours, c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS'))
study_shapes <- c('pSS'=21, 'UC'=22, 'CO'=23, 'TI'=23, 'SLE'=24, 'MS'=25)

### Bar plot of the number of DEGs in each study ###
deg_size_list <- list()
for(i in 1:length(study_list)){
    tmp <- lapply(study_list[[i]], function(x){
        nrow(subset(x, FDR < 0.05 & abs(logFC) > 0.1))
    })
    deg_size_list[[i]] <- data.frame(celltype = replace.names(gsub('_', '.', names(tmp))), size = unlist(tmp))
}
names(deg_size_list) <- names(study_list)

lapply(deg_size_list, function(x){
    sum(x$size)
})

plots_list <- list()
for (name in names(deg_size_list)) {
    df <- deg_size_list[[name]][order(deg_size_list[[name]]$size, decreasing = FALSE),]
    df$celltype <- factor(df$celltype, levels = df$celltype)
    p <- ggplot(df, aes(x = celltype, y = size, fill = celltype)) +
        geom_col() +
        theme_minimal() +
        scale_fill_manual(values = celltypes_colour) +
        # geom_text(aes(label = size), size = 3, hjust = -0.5) +
        xlab('') +
        ylab('# of DEGs') +
        theme(legend.position = 'none') +
        ggtitle(name) +
        coord_flip()
    # Store the plot in the list
    plots_list[[name]] <- p
}

# Save the arranged plots to a single PDF
pdf('Aim_1/combined_deg_size_barplots.pdf', width = 10, height = 10)  # Adjust size as needed
grid.arrange(grobs = plots_list, ncol = 3, nrow = 2)
dev.off()

### Read in DisGeNET ###
pSS_disgene <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/pSS.tsv')$Gene
UC_disgene <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/UC.tsv')$Gene
CO_disgene <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/CD_colon.tsv')$Gene
TI_disgene <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/CD_TI.tsv')$Gene
SLE_disgene <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/SLE.tsv')$Gene
MS_disgene <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/MS.tsv')$Gene

disgene_list <- list('pSS'=pSS_disgene, 'UC'=UC_disgene, 'CO'=CO_disgene, 'TI'=TI_disgene, 'SLE'=SLE_disgene, 'MS'=MS_disgene)

### Calculate chi-squared test for DisGeNET enrichment ###
chisq_list <- list()
for(names in names(study_list)){
    study <- study_list[[names]]
    disgene <- disgene_list[[names]]
    results <- lapply(names(study), function(cell){
        tmp <- chisq.test.edgeR(study[[cell]], disgene, logfc = 0.1)
        size <- nrow(subset(study[[cell]], FDR < 0.05 & abs(logFC) > 0.1 & gene %in% disgene))
        data.frame(celltype = cell, p.value = tmp$p.value, size = size)
    })
    results <- bind_rows(results)
    results$FDR <- p.adjust(results$p.value, method='fdr')
    results$celltype <- replace.names(gsub('_', '.', results$celltype))
    chisq_list[[names]] <- results
}

save(chisq_list, file = 'Aim_1/disgene_enrichment.Rdata')

# Initialize an empty list to store the plots
plots_list <- list()

for (name in names(chisq_list)) {
    p <- ggplot(chisq_list[[name]], aes(x = celltype, y = -log10(FDR), size = size, color = -log10(FDR))) +
        geom_point() +
        scale_color_gradient2(low = 'grey', mid = 'blue', high = 'red') +
        theme_minimal() +
        geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
        xlab('') +
        ylab('-log10(FDR)') +
        ggtitle(name) +
        coord_flip()
    # Store the plot in the list
    plots_list[[name]] <- p
}

# Save the arranged plots to a single PDF
pdf('Aim_1/combined_disgene_enrichment_dotplots.pdf', width = 20, height = 10)  # Adjust size as needed
grid.arrange(grobs = plots_list, ncol = 3, nrow = 2)
dev.off()

### Chi-squared test for enrichment of Immport genes ###    
immport_files <- list.files('/directflow/SCCGGroupShare/projects/lacgra/immport_genelist', full.names=TRUE)
immport <- lapply(immport_files, function(x){
    read.delim(x)$Symbol
})
names(immport) <- gsub('_', ' ', gsub('.txt', '', basename(immport_files)))

### Iterate over each gene set in immport and perform chi-squared test ###
final_list <- list()
for(geneset in names(immport)){
    chisq_list <- list()
    study_names <- c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS')

    for(i in seq_along(study_list)){
        study <- study_list[[i]]
        result <- lapply(names(study), function(cell){
            exp <- study[[cell]]
            if(length(immport[[geneset]] %in% exp$gene) > 0){
                tmp <- chisq.test.edgeR(exp, immport[[geneset]], logfc = 0.1)
                size <- nrow(subset(exp, FDR < 0.05 & abs(logFC) > 0.1 & gene %in% immport[[geneset]]))
                data.frame(celltype = cell, p.value = tmp$p.value, size = size)
            } else {
                data.frame(celltype = cell, p.value = NA, size = 0)
            }
        })
        names(result) <- names(study)
        result.df <- bind_rows(result, .id='celltype')
        result.df$FDR <- p.adjust(result.df$p.value, method='fdr')
        chisq_list[[study_names[i]]] <- result.df
    }
    final_list[[geneset]] <- chisq_list
}

save(final_list, file = 'Aim_1/immport_enrichment.Rdata')

for(geneset in names(final_list)){
    wide_df <- bind_rows(final_list[[geneset]], .id='study')
    mtx <- dcast(wide_df, celltype ~ study, value.var='FDR')
    mtx <- mtx %>%
    mutate_at(vars(-celltype), as.numeric)
    rownames(mtx) <- replace.names(gsub('_', '.', mtx$celltype))
    mtx <- as.matrix(mtx[,-1])
    mtx[is.na(mtx)] <- 1

    pdf(paste0('Aim_1/immport_',geneset, '_heatmap.pdf'))
    col <- colorRamp2(c(0, 5), c('white','red'))
    draw(Heatmap(-log10(mtx), name='FDR', show_row_names=TRUE,
    col = col, show_column_names=TRUE, cluster_columns=TRUE, cluster_rows=TRUE, 
    column_title=geneset, column_title_gp = gpar(fontsize = 12)))
    dev.off()
}

# Plot best gene set on one page
top_pathways <- c("Antigen Processing and Presentation", "Antimicrobials",
"Cytokines", "TCR Signaling Pathway")

plots_list <- list()
for(geneset in top_pathways){
    wide_df <- bind_rows(final_list[[geneset]], .id='study')
    mtx <- dcast(wide_df, celltype ~ study, value.var='FDR')
    mtx <- mtx %>%
    mutate_at(vars(-celltype), as.numeric)
    rownames(mtx) <- replace.names(gsub('_', '.', mtx$celltype))
    mtx <- as.matrix(mtx[,-1])
    mtx[is.na(mtx)] <- 1

    
    col <- colorRamp2(c(0, 5), c('white','red'))
    p <- Heatmap(-log10(mtx), name='FDR', show_row_names=TRUE,
    col = col, show_column_names=TRUE, cluster_columns=TRUE, cluster_rows=TRUE, 
    column_title=geneset, column_title_gp = gpar(fontsize = 12), row_names_gp = gpar(fontsize = 8))
    # Capture the heatmap as a grob
    ht_grob <- grid.grabExpr(draw(p, newpage = FALSE))

    # Store the grob in the list
    plots_list[[geneset]] <- ht_grob
}

pdf('Aim_1/immport_top_pathways_heatmaps.pdf', width = 10, height = 10)
grid.arrange(grobs = plots_list, ncol = 2, nrow = 2)
dev.off()

lapply(top_pathways, function(x){
    tmp <- final_list[[x]]
    tmp2 <- lapply(tmp, function(y){
        subset(y, FDR < 0.05)$celltype
    })
    sort(table(unlist(tmp2)))
})


table(unlist(lapply(study_list, function(x){
    tmp <- subset(x[['DC2']], FDR < 0.05)
    tmp$gene[tmp$gene %in% immport[['Antimicrobials']]]
})))
immport[['Antigen Processing and Presentation']]

### GSEA overlap ###
load('pSS_GSE157278/Aim_1_2024/figure.data/fgsea_list.RData')
pSS_fgsea <- bind_rows(fgsea_list, .id='celltype')
load('UC_GSE182270/Aim_1_2024/figure.data/fgsea_list.RData')
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

fgsea_df$celltype_study <- paste(fgsea_df$celltype, fgsea_df$study, sep=':')

fgsea_wide <- dcast(pathway ~ celltype_study, value.var='NES', data=fgsea_df)
rownames(fgsea_wide) <- gsub('HALLMARK_', '', fgsea_wide$pathway)
fgsea_wide <- fgsea_wide[,-1]
fgsea_wide[is.na(fgsea_wide)] <- 0
study <- gsub('.+:', '', colnames(fgsea_wide))
celltype <- gsub(':.+', '', colnames(fgsea_wide))

pdf('Aim_1/fgsea_heatmap.pdf', width = 10, height = 10)
anno <- HeatmapAnnotation(
  study = study, 
  celltype = celltype, 
  col = list(
    study = c(
      'pSS' = '#E69F00',  # orange
      'UC' = '#009E73',  # light blue
      'CO' = '#56B4E9',  # teal
      'TI' = '#0072B2',  # darker teal
      'SLE' = '#D55E00',  # reddish orange
      'MS' = '#CC79A7'
    ),
    celltype = celltypes_colour[celltype]
  )
)
col <- colorRamp2(c(-2, 0, 2), c('#67a9cf', 'white', '#ef8a62'))  # light blue to white to coral red
Heatmap(
  as.matrix(fgsea_wide), 
  col = col, 
  name = 'NES', 
  show_row_names = TRUE, 
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8), 
  column_names_gp = gpar(fontsize = 5),
  top_annotation = anno
)
dev.off()


### Chi-squared test for enrichment of escapees in DEGs ###
chisq_list <- list()
study_names <- c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS')

for(i in seq_along(study_list)){
    study <- study_list[[i]]
    result <- lapply(names(study), function(cell){
        exp <- study[[cell]]
        if(length(rownames(escape) %in% exp$gene) > 0){
            tmp <- chisq.test.edgeR(exp, rownames(escape), logfc = 0.1)
            size <- nrow(subset(exp, FDR < 0.05 & abs(logFC) > 0.1 & gene %in% rownames(escape)))
            data.frame(celltype = cell, p.value = tmp$p.value, size = size)
        } else {
            data.frame(celltype = cell, p.value = NA, size = 0)
        }
    })
    names(result) <- names(study)
    result.df <- bind_rows(result, .id='celltype')
    result.df$FDR <- p.adjust(result.df$p.value, method='fdr')
    chisq_list[[study_names[i]]] <- result.df
}

save(chisq_list, file = 'Aim_1/escape_enrichment.Rdata')

table(unlist(lapply(chisq_list, function(x){
    subset(x, FDR < 0.05)$celltype
})))

# Initialize an empty list to store the plots
plots_list <- list()

for (name in names(chisq_list)) {
    df <- chisq_list[[name]][order(-log10(chisq_list[[name]]$FDR), decreasing = FALSE),]
    df$celltype <- factor(df$celltype, levels = df$celltype)
    p <- ggplot(df, aes(x = celltype, y = -log10(FDR), size = size, color = -log10(FDR))) +
        geom_point() +
        scale_color_gradient2(low = 'grey', mid = 'blue', high = 'red') +
        theme_minimal() +
        geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
        xlab('') +
        ylab('-log10(FDR)') +
        ggtitle(name) +
        coord_flip()
    # Store the plot in the list
    plots_list[[name]] <- p
}

# Save the arranged plots to a single PDF
pdf('Aim_1/combined_escape_enrichment_dotplots.pdf', width = 20, height = 10)  # Adjust size as needed
grid.arrange(grobs = plots_list, ncol = 3, nrow = 2)
dev.off()

# Find which genes are shared between studies
a <- subset(CO[['Tem_Effector_helper_T_cells']], logFC > 0.1 & FDR < 0.05 & gene %in% rownames(escape))$gene
b <- subset(pSS[['Tem_Effector_helper_T_cells']], logFC > 0.1 & FDR < 0.05 & gene %in% rownames(escape))$gene
intersect(a, b)

a <- subset(CO[['Tem_Effector_helper_T_cells']], logFC < -0.1 & FDR < 0.05 & gene %in% rownames(escape))$gene
b <- subset(pSS[['Tem_Effector_helper_T_cells']], logFC < -0.1 & FDR < 0.05 & gene %in% rownames(escape))$gene
intersect(a, b)

# list escape genes in each study
escape_degs <- lapply(study_list, function(x){
    unique(unlist(lapply(x, function(y){
        subset(y, FDR < 0.05 & abs(logFC) > 0.1 & gene %in% rownames(escape))$gene
    })))
})

upset_mtx <- fromList(escape_degs)
rownames(upset_mtx) <- unique(unlist(escape_degs))

pdf('Aim_1/escape_upset.pdf')
upset(upset_mtx, order.by = 'freq', sets.bar.color = study_colours, nsets = 6)
dev.off()

upset_mtx$sum <- rowSums(upset_mtx)
upset_mtx <- upset_mtx[order(upset_mtx$sum, decreasing = TRUE),]
save(upset_mtx, file = 'Aim_1/escape_upset_mtx.Rdata')


rownames(upset_mtx[1:20,])

upset_mtx[rownames(upset_mtx) %in% X.gwas,]


### Chi-squared test for enrichment of XIST binding proteins ###
### XIST RBP ###
XIST_RBP <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/XIST_RBP.txt')
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

chisq_list <- list()
study_names <- c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS')

for(i in seq_along(study_list)){
    study <- study_list[[i]]
    result <- lapply(names(study), function(cell){
        exp <- study[[cell]]
        if(length(XIST_RBP %in% exp$gene) > 0){
            tmp <- chisq.test.edgeR(exp, XIST_RBP, logfc = 0.1)
            size <- nrow(subset(exp, FDR < 0.05 & abs(logFC) > 0.1 & gene %in% XIST_RBP))
            data.frame(celltype = cell, p.value = tmp$p.value, size = size)
        } else {
            data.frame(celltype = cell, p.value = NA, size = 0)
        }
    })
    names(result) <- names(study)
    result.df <- bind_rows(result, .id='celltype')
    result.df$FDR <- p.adjust(result.df$p.value, method='fdr')
    chisq_list[[study_names[i]]] <- result.df
}

save(chisq_list, file = 'Aim_1/XIST_RBP_enrichment.Rdata')

# Initialize an empty list to store the plots
plots_list <- list()

for (name in names(chisq_list)) {
    df <- chisq_list[[name]][order(-log10(chisq_list[[name]]$FDR), decreasing = FALSE),]
    df$celltype <- factor(df$celltype, levels = df$celltype)
    p <- ggplot(df, aes(x = celltype, y = -log10(FDR), size = size, color = -log10(FDR))) +
        geom_point() +
        scale_color_gradient2(low = 'grey', mid = 'blue', high = 'red') +
        theme_minimal() +
        geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
        xlab('') +
        ylab('-log10(FDR)') +
        ggtitle(name) +
        coord_flip()
    # Store the plot in the list
    plots_list[[name]] <- p
}

# Save the arranged plots to a single PDF
pdf('Aim_1/combined_XIST_RBP_enrichment_dotplots.pdf', width = 20, height = 10)  # Adjust size as needed
grid.arrange(grobs = plots_list, ncol = 3, nrow = 2)
dev.off()

### Overlap of cell types ###
mapply(function(x, y) {
    tmp1 <- subset(x, FDR < 0.05)
    tmp2 <- subset(y, FDR < 0.05)
    intersect(tmp1$celltype, tmp2$celltype)
}, escape_enriched, RBP)

table(unlist(lapply(chisq_list, function(x){
    subset(x, FDR < 0.05)$celltype
})))


lapply(study_list, function(x){
    lapply(x, function(y){
        'XIST' %in% subset(y, FDR < 0.05)$gene
    })
})



### Find common cell types between studies ###
study_celltypes <- lapply(study_list, function(x) names(x))
names(study_celltypes) <- c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS')
celltypes_mtx <- fromList(study_celltypes)
rownames(celltypes_mtx) <- unique(unlist(study_celltypes))
### select common cell types with at least 5 samples ###
common <- rownames(celltypes_mtx)[rowSums(celltypes_mtx) >= 5]


### Identify Jaccard index of upregulated genes between common cell types ###
jaccard_index <- function(set1, set2){
    set1 <- set1[!is.na(set1)]
    set2 <- set2[!is.na(set2)]
    intersection_size <- length(intersect(set1, set2))
    union_size <- length(union(set1, set2))
    return(intersection_size/union_size)
}

jaccard_matrix_list <- list()
for(i in 1:length(common)){
    celltype <- lapply(study_list, function(x) x[[common[i]]])
    degs <- lapply(celltype, function(x){
        if(is.null(x)){
            return(NULL)
        }
        subset(x, abs(logFC) > 0.1 & FDR < 0.05)$gene
        # x$gene
    })

    jaccard_matrix <- matrix(NA, nrow=length(study_list), ncol=length(study_list))
    rownames(jaccard_matrix) <- names(study_list)
    colnames(jaccard_matrix) <- names(study_list)

    for(n in 1:length(study_list)){
        for(m in 1:length(study_list)){
            if(is.null(degs[[n]]) | is.null(degs[[m]])){
                jaccard_matrix[n,m] <- NA
            } else {
            jaccard_matrix[n,m] <- jaccard_index(degs[[n]], degs[[m]])
        }
    }
    jaccard_matrix_list[[common[i]]] <- jaccard_matrix
    }
}

lapply(jaccard_matrix_list, function(x){
    max(x[upper.tri(x, diag=FALSE)], na.rm=TRUE)
})


heatmaps_list <- lapply(names(jaccard_matrix_list), function(x) {
    col <- colorRamp2(c(0, 1), c('blue', 'red'))
    ht <- Heatmap(
        jaccard_matrix_list[[x]],
        name = 'Jaccard index',
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        column_title = replace.names(gsub('_', '.', x)),
        column_title_gp = gpar(fontsize = 5),
        col = col
    )
    return(ht)
})

heatmap_grobs <- lapply(heatmaps_list, function(ht) {
    # Use grid.grabExpr to capture the heatmap as a grob
    grob <- grid.grabExpr(draw(ht, newpage = FALSE))
    return(grob)
})

pdf('Aim_1/jaccard_heatmaps_all_V2.pdf', width=10, height=10)  # Adjust size as needed
grid.arrange(grobs = heatmap_grobs, ncol = 4, nrow = 3)
dev.off()

pdf('Aim_1/jaccard_heatmaps_degs_V2.pdf', width=10, height=10)  # Adjust size as needed
grid.arrange(grobs = heatmap_grobs, ncol = 4, nrow = 3)
dev.off()

### Overlap of GSEA results ###
load('pSS_GSE157278/Aim_1_2024/figure.data/fgsea_list.RData')
pSS_GSEA 


### Find common genes between DEGs and NMF features ###
pSS_NMF_files <- list.files('pSS_GSE157278/NMF', pattern='features.csv', full.names=TRUE)
pSS_NMF <- lapply(pSS_NMF_files, function(x) read.csv(x, header=TRUE)$x)
names(pSS_NMF) <- gsub('_features.csv', '', basename(pSS_NMF_files))
UC_NMF_files <- list.files('UC_GSE182270/NMF', pattern='features.csv', full.names=TRUE)
UC_NMF <- lapply(UC_NMF_files, function(x) read.csv(x, header=TRUE)$x)
names(UC_NMF) <- gsub('_features.csv', '', basename(UC_NMF_files))
CO_NMF_files <- list.files('CD_Kong/colon/NMF', pattern='features.csv', full.names=TRUE)
CO_NMF <- lapply(CO_NMF_files, function(x) read.csv(x, header=TRUE)$x)
names(CO_NMF) <- gsub('_features.csv', '', basename(CO_NMF_files))
TI_NMF_files <- list.files('CD_Kong/TI/NMF', pattern='features.csv', full.names=TRUE)
TI_NMF <- lapply(TI_NMF_files, function(x) read.csv(x, header=TRUE)$x)
names(TI_NMF) <- gsub('_features.csv', '', basename(TI_NMF_files))
SLE_NMF_files <- list.files('lupus_Chun/NMF', pattern='features.csv', full.names=TRUE)
SLE_NMF <- lapply(SLE_NMF_files, function(x) read.csv(x, header=TRUE)$x)
names(SLE_NMF) <- gsub('_features.csv', '', basename(SLE_NMF_files))
MS_NMF_files <- list.files('MS_GSE193770/NMF', pattern='features.csv', full.names=TRUE)
MS_NMF <- lapply(MS_NMF_files, function(x) read.csv(x, header=TRUE)$x)
names(MS_NMF) <- gsub('_features.csv', '', basename(MS_NMF_files))

NMF_list <- list('pSS'=pSS_NMF, 'UC'=UC_NMF, 'CO'=CO_NMF, 'TI'=TI_NMF, 'SLE'=SLE_NMF, 'MS'=MS_NMF)
save(NMF_list, file = 'Aim_1/NMF_list.Rdata')



result <- list()
for(study in names(study_list)){
    degs <- study_list[[study]]
    NMF_features <- NMF_list[[study]]
    tmp <- lapply(names(degs), function(x){
        set1 <- degs[[x]][degs[[x]] %in% chrX]
        set2 <- NMF_features[[x]][NMF_features[[x]] %in% chrX]
        intersect(set1, set2)
    })
    names(tmp) <- names(degs)
    result[[study]] <- tmp
}

foo <- unlist(lapply(result, function(x){
    unlist(x)
}))
plot(density(foo))


jaccard_matrix_list <- list()
for(i in 1:length(common)){
    features <- lapply(celltype, function(x){
        if(is.null(x)){
            return(NULL)
        }
        unlist(x)
    })

    jaccard_matrix <- matrix(NA, nrow=length(NMF_list), ncol=length(NMF_list))
    rownames(jaccard_matrix) <- names(NMF_list)
    colnames(jaccard_matrix) <- names(NMF_list)

    for(n in 1:length(NMF_list)){
        for(m in 1:length(NMF_list)){
            if(is.null(features[[n]]) | is.null(features[[m]])){
                jaccard_matrix[n,m] <- NA
            } else {
                jaccard_matrix[n,m] <- jaccard_index(features[[n]], features[[m]])
            }
        }
    }
    jaccard_matrix_list[[common[i]]] <- jaccard_matrix
}


NMF.files <- list.files('NMF', pattern='features.csv', full.names=TRUE)
NMF.features <- lapply(NMF.files, function(x) read.csv(x, header=FALSE)$V1)
names(NMF.features) <- gsub('_features.csv', '', basename(NMF.files))

lapply(names(degs), function(x){
    length(intersect(degs[[x]], NMF.features[[x]]))/length(degs[[x]]) * 100
})

lapply(names(degs), function(x){
    jaccard_index(degs[[x]], NMF.features[[x]])
})