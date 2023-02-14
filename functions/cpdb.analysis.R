library(ggplot2)
library(reshape2)
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')

result <- read.delim('cpdb/analysis/significant_means.txt')

degs <- edgeR.list('psuedobulk', logfc=0.5)
names(degs) <- gsub('.edgeR-LRT', '', names(degs))

# subset result to only include names(degs)
lapply(strsplit(colnames(result)[-c(1:12)], '\\.'), function(x){
    names(degs) %in% x
})



Fibroblasts <- result[, c(1:12, grep('Fibroblasts', colnames(result)))]
Macrophages <- result[, c(1:12, grep('Macrophages', colnames(result)))]
Plasma_cells <- result[, c(1:12, grep('Plasma_cells', colnames(result)))]
Tem_Trm_cytotoxic_T_cells <- result[, c(1:12, grep('Tem_Trm_cytotoxic_T_cells', colnames(result)))]

cells <- list(Fibroblasts, Macrophages, Plasma_cells, Tem_Trm_cytotoxic_T_cells)

lapply(1:4, function(x){
    df <- cells[[x]]
    res <- df[df$gene_a %in% degs[[x]]$gene | df$gene_b %in% degs[[x]]$gene,]
    apply(res[, 13:ncol(res)], 1, is.na)
})

test <- Macrophages[Macrophages$gene_a %in% degs[[2]]$gene | Macrophages$gene_b %in% degs[[2]]$gene,]
test[,!apply(test, 2, is.na)]