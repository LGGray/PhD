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

### Set colours and shape for celltypes and study
celltypes_colour <- scico(length(common), palette = 'lipari')
celltypes_colour <- setNames(celltypes_colour, replace.names(gsub('_', '.', common)))
study_colours <- scico(6, palette = 'roma')
study_colours <- setNames(study_colours, c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS'))
study_shapes <- c('pSS'=21, 'UC'=22, 'CO'=23, 'TI'=23, 'SLE'=24, 'MS'=25)

deg_size_list <- list()
for(cell in common){
    tmp <- lapply(study_list, function(x){
        if(is.null(x[[cell]])){
            return(0)
        } else {
            return(nrow(subset(x[[cell]], FDR < 0.05 & abs(logFC) > 0.1)))
    }
    })
    deg_size_list[[cell]] <- data.frame(study = names(study_list), size = unlist(tmp))
}
sort(unlist(lapply(deg_size_list, function(x) sum(x$size))))
deg_size <- bind_rows(deg_size_list, .id='celltype')
deg_size$celltype <- replace.names(gsub('_', '.', deg_size$celltype))

pdf('Aim_1/deg_size_barplot.pdf', width=10, height=10)
ggplot(deg_size, aes(x=study, y=size, fill=study)) +
    geom_col() +
    theme_minimal() +
    scale_fill_manual(values=study_colours) +
    geom_text(aes(label=size), vjust=-0.5, size=3) +
    xlab('Study') +
    ylab('# of DEGs') +
    ggtitle('Number of DEGs in common cell types') +
    coord_flip() +
    facet_wrap(~celltype, scales='free_y')
dev.off()

### Chi-squared test for enrichment of X chromosome ###
chisq_list <- list()
study_names <- c('pSS', 'UC', 'CO', 'TI', 'SLE', 'MS')

for(i in seq_along(study_list)){
    study <- study_list[[i]]
    result <- lapply(names(study), function(cell){
        exp <- study[[cell]]
        if(length(chrX %in% exp$gene) > 0){
            tmp <- chisq.test.edgeR(exp, chrX, logfc = 0.1)
            size <- nrow(subset(exp, FDR < 0.05 & abs(logFC) > 0.1 & gene %in% chrX))
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

save(chisq_list, file = 'Aim_1/chrX_enrichment.Rdata')

lapply(chisq_list, function(x){
    nrow(subset(x, FDR < 0.05))
})

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

plot.data <- chisq_list[['SLE']]
plot.data$celltype <- replace.names(gsub('_', '.', plot.data$celltype))
pdf('Aim_1/escape_enrichment_dotplot.pdf')
ggplot(plot.data, aes(x=celltype, y=-log10(FDR), size=size, color=-log10(FDR))) +
    geom_point() +
    scale_size_continuous(range=c(1, 10)) +
    scale_color_gradient2(low='grey', mid='blue', high='red') +
    theme_minimal() +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    xlab('') +
    ylab('-log10(FDR)') +
    ggtitle('SLE: Escapee enrichment in DEGs') +
    coord_flip()
dev.off()

### Chi-squared test for enrichment of XIST binding proteins ###
### XIST RBP ###
XIST_RBP <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/XIST_RBP.txt')
XIST_RBP$Gene.Symbol
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

plot.data <- bind_rows(chisq_list[c('pSS', 'CO', 'TI', 'SLE')], .id='study')
plot.data$celltype <- replace.names(gsub('_', '.', plot.data$celltype))

pdf('Aim_1/XIST_RBP_enrichment_dotplot.pdf', width=10, height=10)
ggplot(plot.data, aes(x=celltype, y=-log10(FDR), size=size, color=-log10(FDR))) +
    geom_point() +
    scale_size_continuous(range=c(1, 10)) +
    scale_color_gradient2(low='grey', mid='blue', high='red') +
    theme_minimal() +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    xlab('') +
    ylab('-log10(FDR)') +
    ggtitle('XIST RBP enrichment in DEGs') +
    coord_flip() +
    facet_wrap(~study, scales='free_y')
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

for(geneset in names(immport)){
    print(geneset)
    print(lapply(final_list[[geneset]], function(x){
    nrow(subset(x, FDR < 0.05))
    }))
}

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
        #x$gene
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


heatmaps <- lapply(names(jaccard_matrix_list), function(x){
    Heatmap(jaccard_matrix_list[[x]], name='Jaccard index', show_row_names=TRUE, show_column_names=TRUE,
    cluster_columns=FALSE, cluster_rows=FALSE, column_title=replace.names(gsub('_', '.', x)), column_title_gp = gpar(fontsize = 5))
})
pdf('Aim_1/jaccard_heatmaps_degs.pdf', width=20, height=5)
heatmaps[[1]] + heatmaps[[2]] + heatmaps[[3]] + heatmaps[[4]] + heatmaps[[5]] + heatmaps[[6]] + heatmaps[[7]] + heatmaps[[8]] + heatmaps[[9]] + heatmaps[[10]] + heatmaps[[11]]
dev.off()

pdf('Aim_1/jaccard_heatmaps_all.pdf', width=20, height=5)
heatmaps[[1]] + heatmaps[[2]] + heatmaps[[3]] + heatmaps[[4]] + heatmaps[[5]] + heatmaps[[6]] + heatmaps[[7]] + heatmaps[[8]] + heatmaps[[9]] + heatmaps[[10]] + heatmaps[[11]]
dev.off()

length(unique(unlist(lapply(MS, function(x) x$gene))))