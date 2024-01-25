library(Seurat)
library(dplyr)
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')


# Script to generate supplementary tables for DEG analysis
pbmc <- readRDS('pbmc.female.RDS')

# Object metrics
df <- pbmc@meta.data %>%
  group_by(cellTypist) %>%
  summarize('N cells' = n(),
            'Samples' = n_distinct(individual),
            'Cases' = sum(grepl('disease', unique(paste(individual, condition)))),
            'Controls' = sum(grepl('control', unique(paste(individual, condition))))) %>%
  as.data.frame()
df$cellTypist <- gsub('/| |-', '_', df$cellTypist)

deg <- deg.list('differential.expression/edgeR', logfc=0.5)
deg.chrX <- lapply(deg, function(x) subset(x, gene %in% rownames(chrX))$gene)

deg.results <- lapply(deg, function(x){
  if(nrow(x) == 0) {
    return(data.frame(Upregulated=0, up.X=0, Downregulated=0, down.X=0))
  } else {
    up.X <- nrow(subset(x, x[,2] > 0.1 & gene %in% rownames(chrX)))
    down.X <- nrow(subset(x, x[,2] < -0.1 & gene %in% rownames(chrX)))
  }
  data.frame(Upregulated=sum(x[,2] > 0.1), 
             up.X=up.X,
             Downregulated=sum(x[,2] < -0.1),
             down.X=down.X)
})

deg.results <- bind_rows(deg.results, .id='celltype')
deg.results$total <- deg.results$Upregulated + deg.results$Downregulated
deg.results$total.X <- deg.results$up.X + deg.results$down.X
deg.results$celltype <- replace.names(gsub('_', '.', deg.results$celltype))
write.table(deg.results, 'differential.expression/deg_cellCount.metrics.csv', sep=',', quote=F, row.names=F)

metrics <- merge(df, deg.results, by='cellTypist')

metrics <- metrics %>%
  mutate(up.X = ifelse(up.X > 0, paste0(up.X, ' (', round(up.X/Upregulated*100, 1), ')'), paste0(up.X, ' (0)')),
       down.X = ifelse(down.X > 0, paste0(down.X, ' (', round(down.X/Downregulated*100, 1), ')'), paste0(down.X, ' (0)')))

colnames(metrics) <- c('Cell type', 'N cells', 'Samples', 'Cases', 'Controls', 
'Upregulated (FDR 0.05, logFC > 0.5)', 'N (%) chrX genes', 
'Downregulated (FDR 0.05, logFC < -0.5)', 'N (%) chrX genes')

write.table(metrics, 'differential.expression/edgeR.metrics.txt', sep='\t', quote=F, row.names=F)

########################################################
features.selected <- function(path, celltype, features){
    tmp <- list.files(path, pattern=paste(celltype, features, 'txt', sep='.'), full.names=T)
    result <- lapply(tmp, function(x) read.delim(x, header=T))
    unique(unlist(result))
}

# if second command line argument is TRUE, generate supplementary tables for ML analysis
if(commandArgs(trailingOnly = TRUE)[2] == 'ML'){
    # Script to generate supplementary tables for ML analysis
    ML <- lapply(names(deg), function(celltype){
        tmp <- readRDS(paste0('exp.matrix/', gsub('_', '.', celltype), '.chrX.RDS'))
        N.cells <- nrow(tmp)
        Samples <- length(unique(tmp$individual))
        Avg.cell.per.sample <- round(N.cells/Samples, 1)
        HVG <- length(features.selected('ML.models/features/', gsub('_', '.', celltype), 'HVG'))
        HVG.chrX <- sum(features.selected('ML.models/features/', gsub('_', '.', celltype), 'HVG') %in% rownames(chrX))
        HVG.DEG <- sum(features.selected('ML.models/features/', gsub('_', '.', celltype), 'HVG') %in% deg[[celltype]]$gene)
        chrX <- length(features.selected('ML.models/features/', gsub('_', '.', celltype), 'chrX'))
        HVG.X <- length(features.selected('ML.models/features/', gsub('_', '.', celltype), 'HVG-X'))
        HVG.random <- length(features.selected('ML.models/features/', gsub('_', '.', celltype), 'HVG-random'))
        return(data.frame(N.cells, Samples, Avg.cell.per.sample, HVG, HVG.chrX, HVG.DEG, chrX, HVG.X, HVG.random))
        })
        names(ML) <- names(deg)

    ML <- bind_rows(ML, .id='Cell type')
    ML <- ML %>% 
    mutate(
        HVG.chrX = ifelse(HVG.chrX == 0, '0 (0)', paste0(HVG.chrX, ' (', round(HVG.chrX/HVG*100, 1), ')')),
        HVG.DEG = ifelse(HVG.DEG == 0, '0 (0)', paste0(HVG.DEG, ' (', round(HVG.DEG/HVG*100, 1), ')'))
    )
    colnames(ML) <- c('Cell type', 'N cells', 'Samples', 
                  'Average cells per samples', 'HVG features selected', 'N (%) features chrX', 
                  'N (%) DE genes overlap', 'chrX features selected', 'HVG-chrX features selected', 'HVG+noise features selected')

    write.table(ML, 'ML.models/selected.features.metrics.txt', sep='\t', quote=F, row.names=F)
} else {
    quit()
}