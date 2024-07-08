library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)

date_threshold <- as.Date("2024-07-04")

methods <- c('boruta', 'enet', 'intersection', 'combined')

metrics_list <- list()
for(i in 1:10){
    metric.files <- unlist(lapply(methods, function(method){
    list.files(paste0('split_', i, '/', method, '/ensemble'), pattern='metrics_', full.names=TRUE)
    }))
    metric.files <- metric.files[as.Date(file.mtime(metric.files)) > date_threshold]
    metrics <- lapply(metric.files, read.csv)
    names(metrics) <- gsub('split_./|ensemble/metrics_|.csv', '', metric.files) %>% gsub('/', '_', .)
    metrics_df <- bind_rows(metrics, .id='celltype') %>%
        mutate(
        gene.set=str_extract(celltype, "HVG.autosome|chrX|autosome|HVG|SLE"),
        method=str_extract(celltype, "boruta|enet|intersection|combined"),
        celltype=gsub('^.+_|.HVG.autosome', '', celltype),
        celltype=gsub('.SLE|.chrX|.autosome|.HVG', '', celltype),
        )
    metrics_list[[i]] <- metrics_df
}
names(metrics_list) <- paste('split', 1:10, sep='_')

### Testing for difference in F1 score between gene sets across splits
compare_gene.set <- lapply(metrics_list, function(x){kruskal.test(F1 ~ gene.set, data=x)$p.value})
compare_gene.set <- t(bind_rows(compare_gene.set, .id='split'))
compare_gene.set <- data.frame(p.value=compare_gene.set[,1], FDR=p.adjust(compare_gene.set[,1], method='fdr'))

### Testing for difference in F1 score between methods across splits
compare_method <- lapply(metrics_list, function(x){kruskal.test(F1 ~ method, data=x)$p.value})
compare_method <- t(bind_rows(compare_method, .id='split'))
compare_method <- data.frame(p.value=compare_method[,1], FDR=p.adjust(compare_method[,1], method='fdr'))

# Combine metrics from all splits
metrics_df <- bind_rows(metrics_list, .id='split')
# Subset for intersection method since no difference in methods
metrics_df <- subset(metrics_df, method == 'intersection')
# Format gene.set as factor
metrics_df$gene.set <- factor(metrics_df$gene.set, levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))
# Format split as factor
metrics_df$split <- factor(metrics_df$split, levels=paste('split', 1:10, sep='_'))

pdf('figures/F1_across_splits.pdf')
ggplot(metrics_df, aes(x=split, y=F1, fill=gene.set)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='split', y='F1 score')
dev.off()

# Average metrics across splits to get a single value for each model
average_metrics <- metrics_df %>% 
    group_by(celltype, gene.set) %>%
    summarise(
        F1=mean(F1),
        n_features=round(mean(n_features),1),
        AUPRC=mean(AUPRC)) %>%
    data.frame()

kruskal.test(F1 ~ gene.set, data=average_metrics)
pairwise.wilcox.test(average_metrics$F1, average_metrics$gene.set, p.adjust.method='BH')

# chrX vs HVG
wilcox.test(average_metrics$F1[average_metrics$gene.set == "chrX"], 
            average_metrics$F1[average_metrics$gene.set == "HVG"], 
            alternative='greater')

# chrX vs HVG.autosome
wilcox.test(average_metrics$F1[average_metrics$gene.set == "chrX"], 
            average_metrics$F1[average_metrics$gene.set == "HVG.autosome"],
            alternative='greater')

# chrX vs SLE
wilcox.test(average_metrics$F1[average_metrics$gene.set == "chrX"], 
            average_metrics$F1[average_metrics$gene.set == "SLE"],
            alternative='greater')

pdf('figures/avg_F1_boxplot.pdf')
ggplot(average_metrics, aes(x=gene.set, y=F1, fill=gene.set)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='Gene Set', y='F1 score')
dev.off()

pdf('figures/avg_AUPRC_boxplot.pdf')
ggplot(average_metrics, aes(x=gene.set, y=AUPRC, fill=gene.set)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='Gene Set', y='AUPRC')
dev.off()

pdf('figures/avg_n_features_boxplot.pdf')
ggplot(average_metrics, aes(x=gene.set, y=n_features, fill=gene.set)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='Gene Set', y='Number of Features')
dev.off()

pdf('F1_by_n_features.pdf')
ggplot(metrics_df, aes(x=log(n_features), y=F1)) +
    geom_point() +
    theme_minimal() +
    labs(x='log(Number of Features)', y='F1 score') +
    geom_smooth(method='lm')
    facet_wrap(~gene.set)
dev.off()


feature_list <- list()
for(i in 1:10){
    feature.files <- list.files(paste0('split_', i, '/features'), pattern='intersection', full.names=TRUE)
    feature.files <- feature.files[as.Date(file.mtime(feature.files)) > date_threshold]
    features <- lapply(feature.files, function(x) read.csv(x)$Feature)
    names(features) <- gsub('^.+_|.csv', '', basename(feature.files))
    feature_list[[i]] <- features
}

celltype.geneset <- names(feature_list[[1]])

result_list <- list()
for(i in celltype.geneset){
    tmp <- lapply(1:10, function(x){
    feature_list[[x]][[i]]
    })
    result_list[[i]] <- tmp
}
names(result_list) <- celltype.geneset

for()

mtx <- fromList(result_list[[1]])
rownames(mtx) <- unique(unlist(result_list[[1]]))
colnames(mtx) <- paste('split', 1:10, sep='_')

col_fun <- colorRamp2(c(0, 1), c("white", "red"))
pdf('figures/test.heatmap.pdf')
Heatmap(as.matrix(mtx), 
        name='Feature', 
        cluster_columns=FALSE, 
        col=col_fun,
        show_heatmap_legend=FALSE)
dev.off()
