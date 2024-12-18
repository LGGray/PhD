library(Seurat)

pbmc <- readRDS('pbmc.female.RDS')

metadata <- unique(pbmc@meta.data[,c('individual', 'condition', 'disease_state')])
rownames(metadata) <- NULL

flare <- subset(metadata, disease_state == 'flare')
colnames(flare) <- c('rownames', 'condition', 'disease_state')

for(i in 1:10){
    test_index <- read.csv(paste0('new_pseudobulk/split_', i, '/test_index.csv'))
    
    test_index$condition <- metadata[match(test_index$rownames, metadata$individual),'condition']
    
    flare_test_index <- subset(test_index, condition == 'control')
    
    flare_test_index <- rbind(flare_test_index, flare[,c('rownames', 'condition')])
    
    write.csv(flare_test_index, paste0('new_pseudobulk/split_', i, '/flare_test_index.csv'), row.names = FALSE)
}

dir.create('new_pseudobulk/flare')

pbmc$development_stage <- as.numeric(gsub('-year-old human stage', '', pbmc$development_stage))

for(cell in levels(pbmc)){
    pbmc.subset <- subset(pbmc, idents = cell)

    # Pseudobulk by average expression
    bulk <- AverageExpression(pbmc.subset, slot='counts', group.by='individual')$RNA
    bulk <- t(bulk)

    # Create metadata to add to pseudobulk matrix
    meta <- unique(pbmc.subset@meta.data[,c('condition', 'individual', 'development_stage', 'ethnicity')])
    meta$cellCount <- sapply(meta$individual, function(id) sum(pbmc.subset$individual == id))
    meta <- meta[match(gsub('_.+', '', rownames(bulk)), meta$individual),]

    bulk <- cbind(class=meta$condition, individual=meta$individual, cellCount=meta$cellCount, age=meta$development_stage,
    ancestry=meta$ethnicity, bulk)
    cell <- gsub('/| |-', '.', cell)
    saveRDS(bulk, paste0('new_pseudobulk/flare/', cell, '.RDS'))
}

### Analysis of flare results ###
library(ggplot2)
library(ggsignif)
library(dplyr)
library(tidyr)
library(stringr)
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')

gene.set.colours <- c('chrX'='#8A0798', 'autosome'='#7DB176', 'HVG'='#44ABD9', 'HVG.autosome'='#1007D9', 'SLE'='#D90750')


metrics_list <- list()
for(i in 1:10){
    metrics.files <- list.files(paste0('split_', i, '/flare'), pattern='metrics', full.names = TRUE)
    metrics <- lapply(metrics.files, read.csv)
    names(metrics) <- gsub('metrics_|.csv', '', basename(metrics.files))
    metrics_df <- bind_rows(metrics, .id='celltype')
    metrics_list[[i]] <- metrics_df
}
metrics <- bind_rows(metrics_list, .id='split')
metrics$gene.set <- str_extract(metrics$celltype, 'chrX|autosome|HVG|HVG.autosome|SLE')
metrics <- metrics %>% 
    mutate(celltype=gsub('^.+_|.HVG.autosome', '', celltype),
        celltype=gsub('.SLE|.chrX|.autosome|.HVG', '', celltype))
metrics$gene.set <- factor(metrics$gene.set, levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))

pdf('figures/flare_boxplot.pdf')
ggplot(metrics, aes(x=gene.set, y=MCC)) +
    geom_boxplot() +
    geom_jitter(width=0.2) +
    # geom_signif(comparisons = list(c('chrX', 'autosome'), c('HVG', 'HVG.autosome'), c('HVG', 'SLE'), c('HVG.autosome', 'SLE')),
    #     map_signif_level=TRUE, textsize=5) +
    labs(x='Gene set', y='MCC')
dev.off()

calc_CI <- function(x){
    mean_score <- mean(x)
    std_dev_score <- sd(x)
    n <- length(x)
    alpha <- 0.05
    # Calculate the t critical value
    t_critical <- qt(alpha / 2, df = n - 1, lower.tail = FALSE)
    # Calculate the margin of error
    margin_error <- t_critical * (std_dev_score / sqrt(n))
    # Calculate the lower and upper bounds of the 95% CI
    lower_bound <- mean_score - margin_error
    upper_bound <- mean_score + margin_error
    return(c(lower_bound, upper_bound))
}

### Average metrics across splits to get a single value for each model ###
average_metrics <- lapply(split(metrics, interaction(metrics$celltype, metrics$gene.set, sep='_')), function(x){
    data.frame(MCC=mean(x$MCC), MCC_lower=calc_CI(x$MCC)[1], MCC_upper=calc_CI(x$MCC)[2],
    n_features=round(mean(x$n_features), 0))
})
average_metrics <- bind_rows(average_metrics, .id='celltype_gene.set') %>%
    separate_wider_delim(celltype_gene.set, delim = "_", names = c('celltype', 'gene.set')) %>%
    data.frame()
average_metrics$gene.set <- factor(average_metrics$gene.set, levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))
average_metrics$celltype <- replace.names(average_metrics$celltype)

kruskal.test(MCC ~ gene.set, data=average_metrics)

metrics$celltype <- replace.names(metrics$celltype)
celltypes <- subset(average_metrics, MCC > 0.7 & gene.set == 'chrX')$celltype
metrics.subset <- subset(metrics, celltype %in% celltypes)

lapply(split(metrics.subset, metrics.subset$celltype), function(x){
    pairwise.wilcox.test(x$MCC, x$gene.set, p.adjust.method='fdr')
})

pdf('figures/avg_flare_MCC_forest.pdf', width=10, height=5)
ggplot(average_metrics, aes(x=MCC, y=celltype)) +
    geom_errorbarh(aes(xmin = MCC_lower, xmax = MCC_upper)) +
    geom_vline(xintercept = 0.7, linetype = 'dotted', color='red') +
    geom_point() +
    # geom_point(data=best_model_point, color='red') +
    labs(x='MCC', y='') +
    theme(axis.text.y=element_text(size=12)) +
    facet_wrap(~gene.set, ncol=5, nrow=1)
dev.off()

pdf('figures/flare_MCC_boxplot.pdf')
comparisons <- list(c('chrX', 'autosome'), c('chrX', 'HVG'), c('chrX', 'SLE'))
y_positions <- c(1.2, 1.4, 1.6)
ggplot(subset(metrics, celltype %in% celltypes), aes(x=gene.set, y=MCC, colour=gene.set)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape=NA, colour='black', fill=NA) +
    geom_signif(
        comparisons=comparisons, 
        map_signif_level=TRUE, 
        test='wilcox.test', 
        color='black', 
        y_position=y_positions # Use dynamic or manual values
    ) +
    theme_minimal() +
    labs(x='', y='MCC') +
    theme(
        axis.text.x=element_blank(),
        strip.text = element_text(size = 10)
    ) +
    scale_colour_manual(values=gene.set.colours, name='Gene Set') +
    facet_wrap(~celltype, strip.position = 'bottom')
dev.off()


