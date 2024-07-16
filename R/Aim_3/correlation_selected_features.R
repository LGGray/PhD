library(dplyr)
library(ggplot2)

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')
gene.set.colours <- c('chrX'='#8A0798', 'autosome'='#7DB176', 'HVG'='#44ABD9', 'HVG.autosome'='#1007D9', 'SLE'='#D90750')

## Testing for correlation between selected features
split_list <- list()
for (split in 1:10){

    feature.files <- list.files(paste0('split_', split, '/features'), pattern='intersection', full.names=TRUE)
    features <- lapply(feature.files, function(x) read.csv(x)$Feature)
    names(features) <- gsub('intersection_|.csv', '', basename(feature.files))

    features <- features[sapply(features, function(x) length(x) > 1)]

    correlation_list <- list()
    for(i in 1:length(features)) {
        mtx <- readRDS(paste0(names(features)[i], '.RDS'))
        mtx <- mtx[, features[[i]]]
        # drop ancestry from mtx
        mtx <- mtx[,colnames(mtx) != 'ancestry']
        correlation <- cor(mtx, method='spearman')
        tmp <- data.frame(
            celltype=names(features)[i],
            min=min(correlation[upper.tri(correlation)]), 
            max=max(correlation[upper.tri(correlation)]), 
            mean=mean(correlation[upper.tri(correlation)])
        )
        correlation_list[[names(features)[i]]] <- tmp
    }
    names(correlation_list) <- names(features)
    correlation_list <- dplyr::bind_rows(correlation_list, .id=NULL)
    split_list[[split]] <- correlation_list
}

correlation_df <- dplyr::bind_rows(split_list, .id='split') %>%
    mutate(
        gene.set=str_extract(celltype, "HVG.autosome|chrX|autosome|HVG|SLE"),
        celltype=gsub('^.+_|.HVG.autosome', '', celltype),
        celltype=gsub('.SLE|.chrX|.autosome|.HVG', '', celltype)
    ) %>%
    data.frame()

correlation_df$gene.set <- factor(correlation_df$gene.set, levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))
correlation_df$celltype <- replace.names(correlation_df$celltype)
write.csv(correlation_df, 'figures/correlation_df.csv')

kruskal.test(max ~ gene.set, data=correlation_df)
pairwise.wilcox.test(correlation_df$max, correlation_df$gene.set, p.adjust.method='fdr')

pdf('figures/feature_correlation.pdf')
ggplot(correlation_df, aes(x=gene.set, y=max, colour=gene.set)) +
    geom_boxplot() +
    theme(axis.text.x=element_blank(),
    strip.text = element_text(size = 4)) +
    labs(x='Gene Set', y='Maximum Correlation') +
    scale_colour_manual(values=gene.set.colours, name='Gene Set') +
    facet_wrap(~celltype, ncol=6, nrow=4)
dev.off()
