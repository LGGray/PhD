library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')

date_threshold <- as.Date("2024-07-04")

gene.set.colours <- c('chrX'='#8A0798', 'autosome'='#7DB176', 'HVG'='#44ABD9', 'HVG.autosome'='#1007D9', 'SLE'='#D90750')
method.colours <- c('boruta'='#0FF5CA', 'enet'='#F59B11', 'intersection'='#1105F2', 'combined'='#F54DEF')
model.colours <- c('logit'='#F4CE03', 'RF'='#F52C2C', 'SVM'='#99B2F5', 'GBM'='#F5B29E', 'MLP'='#26779E')

methods <- c('boruta', 'enet', 'intersection', 'combined')

# Read in metrics for all ML models across gene sets and splits
models_metrics_list <- list()
for(i in 1:10){
    metric.files <- unlist(lapply(methods, function(method){
    list.files(paste0('split_', i, '/', method, '/metrics'), pattern='metrics_', full.names=TRUE)
    }))
    metric.files <- metric.files[as.Date(file.mtime(metric.files)) > date_threshold]
    # Check if metric files are empty
    metric.files <- metric.files[file.size(metric.files) > 0]
    metrics <- lapply(metric.files, read.csv)
    names(metrics) <- gsub('split_./|metrics/|.csv', '', metric.files) %>% 
        gsub('_metrics', '', .) %>%
        gsub('/', '_', .)
    metrics_df <- bind_rows(metrics, .id='celltype') %>%
        mutate(
        model=str_extract(celltype, "logit|RF|SVM|GBM|MLP"),
        gene.set=str_extract(celltype, "HVG.autosome|chrX|autosome|HVG|SLE"),
        method=str_extract(celltype, "boruta|enet|intersection|combined"),
        celltype=gsub('^_.+_', '', celltype),
        celltype=gsub('^.+_|.HVG.autosome', '', celltype),
        celltype=gsub('.SLE|.chrX|.autosome|.HVG', '', celltype),
        )
    models_metrics_list[[i]] <- metrics_df
}
names(models_metrics_list) <- paste('split', 1:10, sep='_')

### Testing for difference in F1 score between method across splits
compare_method <- lapply(models_metrics_list, function(x){kruskal.test(F1 ~ method, data=x)$p.value})
compare_method <- t(bind_rows(compare_method, .id='split'))
compare_method <- data.frame(p.value=compare_method[,1], FDR=p.adjust(compare_method[,1], method='fdr'))

# Combine metrics from all splits
models_metrics_df <- bind_rows(models_metrics_list, .id='split')
models_metrics_df$model <- factor(models_metrics_df$model, levels=c('logit', 'RF', 'SVM', 'GBM', 'MLP'))
models_metrics_df$method <- factor(models_metrics_df$method, levels=c('boruta', 'enet', 'intersection', 'combined'))
models_metrics_df$gene.set <- factor(models_metrics_df$gene.set, levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))
models_metrics_df$split <- factor(models_metrics_df$split, levels=paste('split', 1:10, sep='_'))

pdf('figures/compare_method_across_splits.pdf')
ggplot(models_metrics_df, aes(x=method, y=F1, colour=method, group=method)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='method', y='F1 score') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=method.colours) +
    facet_wrap(~split, ncol=5, nrow=2)
dev.off()

# Filter for intersection method since no difference in methods
models_metrics_list <- lapply(models_metrics_list, function(x) subset(x, method == 'intersection'))

# Testing for difference in F1 score between models across splits
compare_model <- lapply(models_metrics_list, function(x){kruskal.test(F1 ~ model, data=x)$p.value})
compare_model <- t(bind_rows(compare_model, .id='split'))
compare_model <- data.frame(p.value=compare_model[,1], FDR=p.adjust(compare_model[,1], method='fdr'))

# Testing for differences in F1 score between gene sets across splits
compare_gene.set <- lapply(models_metrics_list, function(x){kruskal.test(F1 ~ gene.set, data=x)$p.value})
compare_gene.set <- t(bind_rows(compare_gene.set, .id='split'))
compare_gene.set <- data.frame(p.value=compare_gene.set[,1], FDR=p.adjust(compare_gene.set[,1], method='fdr'))

lapply(models_metrics_list, function(x){
    pairwise.wilcox.test(x$F1, x$gene.set, p.adjust.method='fdr', alternative='greater')
})

pdf('figures/compare_model_across_splits.pdf')
ggplot(models_metrics_df, aes(x=model, y=F1, colour=model, group=model)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='model', y='F1 score') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=model.colours) +
    facet_wrap(~split, ncol=5, nrow=2)
dev.off()

pdf('figures/compare_method_across_splits.pdf')
ggplot(models_metrics_df, aes(x=method, y=F1, colour=method, group=method)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='method', y='F1 score') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=method.colours) +
    facet_wrap(~split, ncol=5, nrow=2)
dev.off()

pdf('figures/compare_geneset_across_splits.pdf')
ggplot(models_metrics_df, aes(x=gene.set, y=F1, colour=gene.set, group=gene.set)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='gene set', y='F1 score') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours) +
    facet_wrap(~split, ncol=5, nrow=2)
dev.off()

# Average metrics across splits to get a single value for each model
average_model_metrics_df <- models_metrics_df %>%
    group_by(celltype, gene.set, model) %>%
    summarise(
        F1=mean(F1),
        n_features=round(mean(n_features),1),
        AUPRC=mean(AUPRC)) %>%
    data.frame()
average_model_metrics_df$model <- factor(average_model_metrics_df$model, levels=c('logit', 'RF', 'SVM', 'GBM', 'MLP'))
average_model_metrics_df$gene.set <- factor(average_model_metrics_df$gene.set, levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))

pdf('figures/avg_model_F1_geneset.pdf')
ggplot(average_model_metrics_df, aes(x=gene.set, y=F1, colour=gene.set)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='Gene Set', y='F1 score') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours)
dev.off()

kruskal.test(F1 ~ gene.set, data=average_model_metrics_df)
pairwise.wilcox.test(average_model_metrics_df$F1, average_model_metrics_df$gene.set, 
    p.adjust.method='fdr')

pdf('figures/avg_model_F1_model.pdf')
ggplot(average_model_metrics_df, aes(x=model, y=F1, colour=model)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='Model', y='F1 score') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=model.colours)
dev.off()

kruskal.test(F1 ~ model, data=average_model_metrics_df)
pairwise.wilcox.test(average_model_metrics_df$F1, average_model_metrics_df$model, 
    p.adjust.method='fdr')

pdf('figures/compare_average_model.pdf')
ggplot(average_model_metrics_df, aes(x=gene.set, y=F1, fill=gene.set)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='gene set', y='F1 score') +
    theme(axis.text.x=element_blank()) +
    scale_fill_manual(values=gene.set.colours) +
    facet_wrap(~model, ncol=5, nrow=1)
dev.off()

### Read in metrics file for ensemble models across gene sets and splits ###
ensmble_metrics_list <- list()
for(i in 1:10){
    metric.files <- list.files(paste0('split_', i, '/intersection/ensemble'), pattern='metrics_', full.names=TRUE)
    metric.files <- metric.files[as.Date(file.mtime(metric.files)) > date_threshold]
    metrics <- lapply(metric.files, read.csv)
    names(metrics) <- gsub('split_./|ensemble/metrics_|.csv', '', metric.files) %>% gsub('/', '_', .)
    metrics_df <- bind_rows(metrics, .id='celltype') %>%
        mutate(
        gene.set=str_extract(celltype, "HVG.autosome|chrX|autosome|HVG|SLE"),
        celltype=gsub('^.+_|.HVG.autosome', '', celltype),
        celltype=gsub('.SLE|.chrX|.autosome|.HVG', '', celltype))
    ensmble_metrics_list[[i]] <- metrics_df
}
names(ensmble_metrics_list) <- paste('split', 1:10, sep='_')

### Testing for difference in F1 score between gene sets across splits
compare_gene.set <- lapply(ensmble_metrics_list, function(x){kruskal.test(F1 ~ gene.set, data=x)$p.value})
compare_gene.set <- t(bind_rows(compare_gene.set, .id='split'))
compare_gene.set <- data.frame(p.value=compare_gene.set[,1], FDR=p.adjust(compare_gene.set[,1], method='fdr'))

# Combine metrics from all splits
ensemble_metrics_df <- bind_rows(ensmble_metrics_list, .id='split')
# Format gene.set as factor
ensemble_metrics_df$gene.set <- factor(ensemble_metrics_df$gene.set, levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))
# Format split as factor
ensemble_metrics_df$split <- factor(ensemble_metrics_df$split, levels=paste('split', 1:10, sep='_'))
# Format celltype names
ensemble_metrics_df$celltype <- replace.names(ensemble_metrics_df$celltype)

pdf('figures/ensemble_geneset_across_splits.pdf')
ggplot(ensemble_metrics_df, aes(x=gene.set, y=F1, colour=gene.set)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='gene set', y='F1 score') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours) +
    facet_wrap(~split, ncol=5, nrow=2)
dev.off()

compare_celltype <- lapply(split(ensemble_metrics_df, ensemble_metrics_df$celltype), function(x){
    kruskal.test(F1 ~ gene.set, data=x)$p.value
})
compare_celltype <- t(bind_rows(compare_celltype, .id='celltype'))
compare_celltype <- data.frame(p.value=compare_celltype[,1], FDR=p.adjust(compare_celltype[,1], method='fdr'))

pdf('figures/ensemble_geneset_celltype.pdf')
ggplot(ensemble_metrics_df, aes(x=gene.set, y=F1, colour=gene.set)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='gene set', y='F1 score') +
    theme(
        axis.text.x=element_blank(),
        strip.text = element_text(size = 4)) +
    scale_colour_manual(values=gene.set.colours) +
    facet_wrap(~celltype, ncol=6, nrow=4)
dev.off()

nonsig_celltypes <- rownames(compare_celltype[compare_celltype$FDR > 0.05,])
lapply(split(ensemble_metrics_df, ensemble_metrics_df$celltype)[nonsig_celltypes], function(x)
    pairwise.wilcox.test(x$F1, x$gene.set, p.adjust.method='fdr')
)

# Identify which split had the highest F1 score for FACS celltypes
ensemble_metrics_df %>%
    filter(gene.set == 'chrX') %>%
    filter(celltype %in% c("Memory.B.cells", "Non.classical.monocytes", 
    "Tcm.Naive.helper.T.cells", "Regulatory.T.cells")) %>%
    group_by(celltype) %>%
    # select the top model based on F1 score
    sort(celltype) %>%
    data.frame()

calc_CI <- function(x){
    mean_F1 <- mean(x)
    std_dev_F1 <- sd(x)
    n <- length(x)
    alpha <- 0.05
    # Calculate the t critical value
    t_critical <- qt(alpha / 2, df = n - 1, lower.tail = FALSE)
    # Calculate the margin of error
    margin_error <- t_critical * (std_dev_F1 / sqrt(n))
    # Calculate the lower and upper bounds of the 95% CI
    lower_bound <- mean_F1 - margin_error
    upper_bound <- mean_F1 + margin_error
    return(c(lower_bound, upper_bound))
}

### Average metrics across splits to get a single value for each model ###
average_ensemble_metrics <- lapply(split(ensemble_metrics_df, interaction(ensemble_metrics_df$celltype, ensemble_metrics_df$gene.set, sep='_')), function(x){
    data.frame(F1=mean(x$F1), F1_lower=calc_CI(x$F1)[1], F1_upper=calc_CI(x$F1)[2],
    AUPRC=mean(x$AUPRC), n_features=round(mean(x$n_features), 0))
})
average_ensemble_metrics <- bind_rows(average_ensemble_metrics, .id='celltype_gene.set') %>%
    separate_wider_delim(celltype_gene.set, delim = "_", names = c('celltype', 'gene.set')) %>%
    data.frame()
average_ensemble_metrics$gene.set <- factor(average_ensemble_metrics$gene.set, levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))

write.csv(average_ensemble_metrics, 'figures/average_ensemble_metrics.csv')

pdf('figures/avg_ensemble_F1_forest.pdf', width=10, height=5)
ggplot(average_ensemble_metrics, aes(x=F1, y=celltype)) +
    geom_point() +
    geom_errorbarh(aes(xmin = F1_lower, xmax = F1_upper)) +
    geom_vline(xintercept = 0.8, linetype = 'dotted', color='red') +
    labs(x='F1-score', y='') +
    facet_wrap(~gene.set, ncol=5, nrow=1)
dev.off()

pdf('figures/avg_ensemble_F1_by_n_features.pdf', width=10, height=5)
ggplot(average_ensemble_metrics, aes(x=log(n_features), y=F1)) +
    geom_point() +
    theme_minimal() +
    labs(x='log(Number of Features)', y='F1 score') +
    geom_smooth(method='lm') +
    stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "pearson") +
    facet_wrap(~gene.set, ncol=5, nrow=1)
dev.off()

pdf('figures/avg_ensemble_n_features_geneset.pdf')
ggplot(average_ensemble_metrics, aes(x=gene.set, y=n_features, colour=gene.set)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='Gene Set', y='Number of Features') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours)
dev.off()

pairwise.wilcox.test(average_ensemble_metrics$n_features, average_ensemble_metrics$gene.set, 
p.adjust.method='fdr')

# Read in the features selected across splits to visualise concordance
feature_list <- list()
for(i in 1:10){
    feature.files <- list.files(paste0('split_', i, '/features'), pattern='intersection', full.names=TRUE)
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

if(!dir.exists('figures/feature_heatmap') == TRUE {dir.create('figures/feature_heatmap')})

for(file in celltype.geneset){

    mtx <- fromList(result_list[[file]])
    rownames(mtx) <- unique(unlist(result_list[[file]]))
    colnames(mtx) <- paste('split', 1:10, sep='_')
    # Order rows by total number of features selected
    mtx <- mtx[order(rowSums(mtx), decreasing = TRUE),]

    col_fun <- colorRamp2(c(0, 1), c("white", "black"))
    pdf(paste0('figures/feature_heatmap/', file, '.pdf'))
    p <- Heatmap(as.matrix(mtx), 
            name='Features', 
            cluster_columns=FALSE,
            cluster_rows=FALSE,
            col=col_fun,
            row_names_gp = gpar(fontsize = 5),
            show_heatmap_legend=FALSE)
    print(p)
    dev.off()
}

### calculate jaccard index between gene sets ###
# Jaccard index
jaccard_index <- function(x, y){
    intersect <- length(intersect(x, y))
    union <- length(union(x, y))
    return(intersect / union)
}

# Pull out features selected across any of the splits
top_features <- list()
for(file in celltype.geneset){
    features <- unique(unlist(result_list[[file]]))
    top_features[[file]] <- features
}

celltypes <- gsub('.HVG.autosome', '', celltype.geneset) %>%
    gsub('.SLE|.chrX|.autosome|.HVG', '', .) %>%
    unique()

jaccard_list <- list()
for( cell in celltypes){
    chrX <- top_features[[paste0(cell, '.chrX')]]
    autosome <- top_features[[paste0(cell, '.autosome')]]
    HVG <- top_features[[paste0(cell, '.HVG')]]
    HVG.autosome <- top_features[[paste0(cell, '.HVG.autosome')]]
    SLE <- top_features[[paste0(cell, '.SLE')]]

    combinations <- combn(c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'), 2)
    results <- data.frame()
    for (i in 1:ncol(combinations)) {
        jaccard_value <- jaccard_index(get(combinations[,i][1]), get(combinations[,][2]))
        results <- rbind(results, data.frame(celltype = cell, set1 = combinations[,i][1],
        set2 = combinations[,i][2], jaccard_index = jaccard_value))
    }

    jaccard_list[[cell]] <- results
}

dir.create('figures/jaccard_index')

pdf('figures/jaccard_index/Age.associated.B.cells.pdf')
ggplot(jaccard_list[['Age.associated.B.cells']], aes(x=set1, y=set2, fill=jaccard_index)) +
    geom_tile() +
    scale_fill_gradient(low='white', high='red') +
    theme_minimal() +
    labs(x='', y='', fill='Jaccard Index') +
    theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

jaccard_df <- bind_rows(jaccard_list, .id=NULL)