library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)
library(ggsignif)
library(enrichR)

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')

gene.set.colours <- c('chrX'='#8A0798', 'autosome'='#7DB176', 'HVG'='#44ABD9', 'HVG.autosome'='#1007D9', 'SLE'='#D90750')
method.colours <- c('boruta'='#BFFB00', 'enet'='#B875B1', 'intersection'='#D2EDF6', 'combined'='#4DB748')
model.colours <- c('logit'='#F4CE03', 'RF'='#BCEA9D', 'SVM'='#99B2F5', 'GBM'='#F5B29E', 'MLP'='#26779E', 'ensemble'='#F5A2F5')

methods <- c('boruta', 'enet', 'intersection', 'combined')
models <- c('logit', 'RF', 'SVM', 'GBM', 'MLP')

# Read in metrics for all ML models across gene sets and splits
models_metrics_list <- list()
for(i in 1:10){
    metric.files <- unlist(lapply(methods, function(method){
    list.files(paste0('split_', i, '/', method, '/metrics'), pattern='metrics_', full.names=TRUE)
    }))
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

# Calculate Mathews Correlation Coefficient (MCC) for each model
MCC <- function(mtx) {
    TP <- mtx[2, 2]
    TN <- mtx[1, 1]
    FP <- mtx[1, 2]
    FN <- mtx[2, 1]
    return((TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
}

methods <- c('boruta', 'enet', 'intersection', 'combined')
for(split in 1:10){
    files <- unlist(lapply(methods, function(method){
    list.files(paste0('split_', split, '/', method, '/metrics'), pattern='confusion.+csv', full.names=TRUE)
    }))
    confusion_list <- lapply(files, read.csv)
    names(confusion_list) <- gsub('split_./|metrics/|_confusion|.csv', '', files)
    result <- lapply(confusion_list, MCC)
    df <- data.frame(celltype=names(result), MCC=unlist(result))
    df <- df %>% mutate(
        model=str_extract(celltype, "logit|RF|SVM|GBM|MLP"),
        gene.set=str_extract(celltype, "HVG.autosome|chrX|autosome|HVG|SLE"),
        method=str_extract(celltype, "boruta|enet|intersection|combined"),
        celltype=gsub('^_.+_', '', celltype),
        celltype=gsub('^.+_|.HVG.autosome', '', celltype),
        celltype=gsub('.SLE|.chrX|.autosome|.HVG', '', celltype),
        )
    df$MCC[!is.finite(df$MCC)] <- 0
    tmp <- models_metrics_list[[split]]
    # merge tmp with MCC
    tmp <- tmp %>% right_join(df, by=c('model'='model', 'gene.set'='gene.set', 'method'='method', 'celltype'='celltype'))
    models_metrics_list[[split]] <- tmp
}

### Testing for difference in MCC score between method across splits
compare_method <- lapply(models_metrics_list, function(x){kruskal.test(MCC ~ method, data=x)$p.value})
compare_method <- t(bind_rows(compare_method, .id='split'))
compare_method <- data.frame(p.value=compare_method[,1], FDR=p.adjust(compare_method[,1], method='fdr'))

# Combine metrics from all splits
models_metrics_df <- bind_rows(models_metrics_list, .id='split')
models_metrics_df$model <- factor(models_metrics_df$model, levels=c('logit', 'RF', 'SVM', 'GBM', 'MLP'))
models_metrics_df$method <- factor(models_metrics_df$method, levels=c('boruta', 'enet', 'intersection', 'combined'))
models_metrics_df$gene.set <- factor(models_metrics_df$gene.set, levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))
models_metrics_df$split <- factor(models_metrics_df$split, levels=paste('split', 1:10, sep='_'))

save(models_metrics_df, file='figures/models_metrics_df.RData')

pdf('figures/compare_method_across_splits.pdf')
ggplot(models_metrics_df, aes(x=method, y=MCC, colour=method, group=method)) +
    geom_jitter(width=0.2, alpha=1) +
    geom_boxplot(outlier.shape = NA, color = 'black', fill=NA) +
    theme_minimal() +
    labs(x='Method', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=method.colours, name='Method') +
    facet_wrap(~split, ncol=5, nrow=2)
dev.off()

# Filter for intersection method since no difference in methods
models_metrics_list <- lapply(models_metrics_list, function(x) subset(x, method == 'intersection'))

# Testing for differences in MCC score between gene sets across splits
compare_gene.set <- lapply(models_metrics_list, function(x){kruskal.test(MCC ~ gene.set, data=x)$p.value})
compare_gene.set <- t(bind_rows(compare_gene.set, .id='split'))
compare_gene.set <- data.frame(p.value=compare_gene.set[,1], FDR=p.adjust(compare_gene.set[,1], method='fdr'))

pdf('figures/compare_geneset_across_splits.pdf')
ggplot(models_metrics_df, aes(x=gene.set, y=MCC, colour=gene.set, group=gene.set)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='Gene Set', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours, name='Gene Set') +
    facet_wrap(~split, ncol=5, nrow=2)
dev.off()

# Testing for difference in MCC score between models across splits
compare_model <- lapply(models_metrics_list, function(x){kruskal.test(MCC ~ model, data=x)$p.value})
compare_model <- t(bind_rows(compare_model, .id='split'))
compare_model <- data.frame(p.value=compare_model[,1], FDR=p.adjust(compare_model[,1], method='fdr'))

lapply(models_metrics_list, function(x){
    pairwise.wilcox.test(x$F1, x$gene.set, p.adjust.method='fdr', alternative='greater')
})

pdf('figures/compare_model_across_splits.pdf')
ggplot(models_metrics_df, aes(x=model, y=MCC, colour=model, group=model)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='Model', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=model.colours, name='Model') +
    facet_wrap(~split, ncol=5, nrow=2)
dev.off()

# Average metrics across splits to get a single value for each model
average_model_metrics_df <- models_metrics_df %>%
    group_by(celltype, gene.set, model) %>%
    summarise(
        MCC=mean(MCC),
        n_features=round(mean(n_features),0)) %>%
    data.frame()
average_model_metrics_df$model <- factor(average_model_metrics_df$model, levels=c('logit', 'RF', 'SVM', 'GBM', 'MLP'))
average_model_metrics_df$gene.set <- factor(average_model_metrics_df$gene.set, levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))

write.csv(average_model_metrics_df,'figures/average_model_metrics.csv')

comparisons <- list(c('chrX', 'autosome'), c('chrX', 'HVG'), c('chrX', 'HVG.autosome'), c('chrX', 'SLE'))
pdf('figures/avg_model_MCC_geneset.pdf')
ggplot(average_model_metrics_df, aes(x=gene.set, y=MCC, colour=gene.set)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape=NA, color='black', fill=NA) +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, test='wilcox.test', color='black') +
    theme_minimal() +
    labs(x='Gene Set', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours, name='Gene Set')
dev.off()

kruskal.test(MCC ~ model, data=average_model_metrics_df)

kruskal.test(MCC ~ model, data=average_model_metrics_df)
MCC_model <- pairwise.wilcox.test(average_model_metrics_df$MCC, average_model_metrics_df$model, 
    p.adjust.method='fdr')
MCC_geneset$p.value[is.na(MCC_geneset$p.value)] <- 1

pdf('figures/pairwise_wilcox_geneset.pdf')
Heatmap(F1_geneset$p.value, col = circlize::colorRamp2(c(0.05, 0.001), c("white", "red")), name = 'FDR')
dev.off()

kruskal.test(MCC ~ model, data=average_model_metrics_df)
pairwise.wilcox.test(average_model_metrics_df$MCC, average_model_metrics_df$model, 
    p.adjust.method='fdr')

comparisons <- list(c('logit', 'MLP'), c('RF', 'MLP'), c('SVM', 'MLP'), c('GBM', 'MLP'))
pdf('figures/avg_model_MCC_model.pdf')
ggplot(average_model_metrics_df, aes(x=model, y=MCC, colour=model)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, test='wilcox.test', color='black') +
    theme_minimal() +
    labs(x='Model', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=model.colours, name='Model')
dev.off()

pdf('figures/compare_average_model.pdf')
ggplot(average_model_metrics_df, aes(x=gene.set, y=AUPRC, fill=gene.set)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='gene set', y='AUPRC') +
    theme(axis.text.x=element_blank()) +
    scale_fill_manual(values=gene.set.colours) +
    facet_wrap(~model, ncol=5, nrow=1)
dev.off()

### Read in metrics file for ensemble models across gene sets and splits ###
ensmble_metrics_list <- list()
for(i in 1:10){
    metric.files <- list.files(paste0('split_', i, '/intersection/ensemble'), pattern='metrics_', full.names=TRUE)
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

# Add MCC to ensemble metrics
for(split in 1:10){
    files <- list.files(paste0('split_', split, '/intersection/ensemble'), pattern='confusion.+csv', full.names=TRUE)
    confusion_list <- lapply(files, read.csv)
    names(confusion_list) <- gsub('confusion_|.csv', '', basename(files))
    result <- lapply(confusion_list, MCC)
    df <- data.frame(celltype=names(result), MCC=unlist(result))
    df <- df %>% mutate(
        gene.set=str_extract(celltype, "HVG.autosome|chrX|autosome|HVG|SLE"),
        celltype=gsub('.HVG.autosome', '', celltype),
        celltype=gsub('.SLE|.chrX|.autosome|.HVG', '', celltype))
    df$MCC[!is.finite(df$MCC)] <- 0
    tmp <- ensmble_metrics_list[[split]]
    # merge tmp with MCC
    tmp <- tmp %>% right_join(df, by=c('gene.set'='gene.set', 'celltype'='celltype'))
    ensmble_metrics_list[[split]] <- tmp
}

### Testing for difference in MCC score between gene sets across splits
compare_gene.set <- lapply(ensmble_metrics_list, function(x){kruskal.test(MCC ~ gene.set, data=x)$p.value})
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

write.csv(ensemble_metrics_df, 'figures/ensemble_metrics.csv', row.names=FALSE)

comparisons <- list(c('chrX', 'autosome'), c('chrX', 'HVG'), c('chrX', 'HVG.autosome'), c('chrX', 'SLE'))
pdf('figures/ensemble_geneset_across_splits_MCC.pdf')
ggplot(ensemble_metrics_df, aes(x=gene.set, y=MCC, colour=gene.set)) +
    geom_boxplot() +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, test='wilcox.test', aes(color='black')) +
    theme_minimal() +
    labs(x='gene set', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours) +
    facet_wrap(~split, ncol=5, nrow=2)
dev.off()

kruskal.test(MCC ~ gene.set, data=ensemble_metrics_df)
pairwise.wilcox.test(ensemble_metrics_df$MCC, ensemble_metrics_df$gene.set, p.adjust.method='fdr')

pdf('figures/ensemble_geneset_MCC.pdf')
comparisons <- list(c('chrX', 'autosome'), c('chrX', 'HVG'), c('chrX', 'HVG.autosome'), c('chrX', 'SLE'))
ggplot(ensemble_metrics_df, aes(x=gene.set, y=MCC, colour=gene.set)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape=NA, colour='black', fill=NA) +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, test='wilcox.test', color='black') +
    theme_minimal() +
    labs(x='Gene Set', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours, name='Gene Set')
dev.off()

compare_celltype <- lapply(split(ensemble_metrics_df, ensemble_metrics_df$celltype), function(x){
    kruskal.test(MCC ~ gene.set, data=x)$p.value
})
compare_celltype <- t(bind_rows(compare_celltype, .id='celltype'))
compare_celltype <- data.frame(p.value=compare_celltype[,1], FDR=p.adjust(compare_celltype[,1], method='fdr'))
subset(compare_celltype, FDR > 0.05)

comparisons <- list(c('chrX', 'autosome'), c('chrX', 'HVG'), c('chrX', 'HVG.autosome'), c('chrX', 'SLE'))
pdf('figures/ensemble_geneset_celltype_MCC.pdf', width=10, height=10)
ggplot(ensemble_metrics_df, aes(x=gene.set, y=MCC, colour=gene.set)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape=NA, colour='black', fill=NA) +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, 
    test='wilcox.test', color='black') +
    theme_minimal() +
    labs(x='Gene Set', y='MCC') +
    theme(
        axis.text.x=element_blank(),
        strip.text = element_text(size = 8)) +
    scale_colour_manual(values=gene.set.colours, name='Gene Set') +
    facet_wrap(~celltype, ncol=5, nrow=5, strip.position = 'bottom')
dev.off()

nonsig_celltypes <- rownames(compare_celltype[compare_celltype$FDR > 0.05,])
lapply(split(ensemble_metrics_df, ensemble_metrics_df$celltype), function(x)
    pairwise.wilcox.test(x$MCC, x$gene.set, p.adjust.method='fdr'))


pdf('figures/top_models_geneset.pdf')
top_models <- subset(ensemble_metrics_df, celltype %in% c('Non-classical monocytes', 'Memory B cells', 'Tcm/Naive helper T cells', 'Regulatory T cells'))
ggplot(top_models, aes(x=gene.set, y=F1, colour=gene.set)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape = NA, color='black', fill=NA) +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, test='wilcox.test', color='black') +
    theme_minimal() +
    labs(x='', y='F1 score') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours, name='Gene Set') +
    facet_wrap(~celltype, strip.position = 'bottom', ncol=2, nrow=2)
dev.off()

# Identify which split had the highest F1 score for FACS celltypes
ensemble_metrics_df %>%
    filter(gene.set == 'chrX') %>%
    filter(celltype %in% c("Memory B cells", "Non classical monocytes", 
    "Tcm Naive helper T cells", "Regulatory T cells")) %>%
    group_by(celltype) %>%
    # select the top model based on F1 score
    sort(celltype) %>%
    data.frame()


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
average_ensemble_metrics <- lapply(split(ensemble_metrics_df, interaction(ensemble_metrics_df$celltype, ensemble_metrics_df$gene.set, sep='_')), function(x){
    data.frame(F1=mean(x$F1), F1_lower=calc_CI(x$F1)[1], F1_upper=calc_CI(x$F1)[2],
    AUPRC=mean(x$AUPRC), AUPRC_lower=calc_CI(x$AUPRC)[1], AUPRC_upper=calc_CI(x$F1)[2], 
    MCC=mean(x$MCC), MCC_lower=calc_CI(x$MCC)[1], MCC_upper=calc_CI(x$MCC)[2],
    n_features=round(mean(x$n_features), 0))
})
average_ensemble_metrics <- bind_rows(average_ensemble_metrics, .id='celltype_gene.set') %>%
    separate_wider_delim(celltype_gene.set, delim = "_", names = c('celltype', 'gene.set')) %>%
    data.frame()
average_ensemble_metrics$gene.set <- factor(average_ensemble_metrics$gene.set, levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))

write.csv(average_ensemble_metrics, 'figures/average_ensemble_metrics.csv', row.names=FALSE)

pdf('figures/avg_ensemble_MCC_geneset.pdf')
comparisons <- list(c('chrX', 'autosome'), c('chrX', 'HVG'), c('chrX', 'HVG.autosome'), c('chrX', 'SLE'))
ggplot(average_ensemble_metrics, aes(x=gene.set, y=MCC, colour=gene.set)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape = NA, colour='black', fill=NA) +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, test='wilcox.test', color='black') +
    theme_minimal() +
    labs(x='Gene Set', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours, name='Gene Set')
dev.off()

pairwise.wilcox.test(average_ensemble_metrics$F1, average_ensemble_metrics$gene.set, p.adjust.method='fdr')

pdf('figures/avg_ensemble_MCC_forest.pdf', width=10, height=5)
ggplot(average_ensemble_metrics, aes(x=MCC, y=celltype)) +
    geom_errorbarh(aes(xmin = MCC_lower, xmax = MCC_upper)) +
    geom_vline(xintercept = 0.7, linetype = 'dotted', color='red') +
    geom_point() +
    # geom_point(data=best_model_point, color='red') +
    labs(x='MCC', y='') +
    theme(axis.text.y=element_text(size=12)) +
    facet_wrap(~gene.set, ncol=5, nrow=1)
dev.off()

subset(average_ensemble_metrics, gene.set == 'chrX' & MCC > 0.7)[,c('celltype', 'MCC')]

pdf('figures/avg_ensemble_MCC_by_n_features.pdf', width=10, height=5)
ggplot(average_ensemble_metrics, aes(x=n_features, y=MCC)) +
    geom_point() +
    theme_minimal() +
    labs(x='Number of Features', y='MCC') +
    geom_smooth(method='lm', se=FALSE) +
    stat_cor(r.digits = 2, cor.coef.name='rho', method = "spearman") +
    facet_wrap(~gene.set, ncol=5, nrow=1, scales='free_x')
dev.off()

comparisons <- list(c('chrX', 'SLE'), c('autosome', 'SLE'), c('HVG', 'SLE'), c('HVG.autosome', 'SLE'))
pdf('figures/avg_ensemble_n_features_geneset.pdf')
ggplot(average_ensemble_metrics, aes(x=gene.set, y=n_features, colour=gene.set)) +
    geom_boxplot() +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, test='wilcox.test', aes(color='black')) +
    theme_minimal() +
    labs(x='Gene Set', y='Number of Features') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours, name='Gene Set')
dev.off()

pairwise.wilcox.test(average_ensemble_metrics$n_features, average_ensemble_metrics$gene.set, 
p.adjust.method='fdr')

# Add everage models and average ensemble to see if ensemble is better
average_ensemble_metrics$model <- 'ensemble'
combined_metrics <- bind_rows(average_model_metrics_df, average_ensemble_metrics)
combined_metrics$model <- factor(combined_metrics$model, levels=c('logit', 'RF', 'SVM', 'GBM', 'MLP', 'ensemble'))

pairwise.wilcox.test(combined_metrics$F1, combined_metrics$model, p.adjust.method='fdr')

pdf('figures/models_ensemble_boxplot.pdf')
combinations <- list(c('logit', 'MLP'), c('RF', 'MLP'), c('SVM', 'MLP'), c('GBM', 'MLP'), c('ensemble', 'MLP'))
ggplot(combined_metrics, aes(x=model, y=F1, colour=model)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape = NA, colour = 'black', fill = NA) +
    geom_signif(comparisons=combinations, map_signif_level=TRUE, test='wilcox.test', color='black') +
    theme_minimal() +
    labs(x='Model', y='F1 score') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=model.colours, name='Model')
dev.off()

kruskal.test(F1 ~ model, data=combined_metrics)
pairwise.wilcox.test(combined_metrics$F1, combined_metrics$model, p.adjust.method='fdr')

### Read in the features selected across splits to visualise concordance ###
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

all_features <- lapply(result_list, function(x) unique(unlist(x)))
top_features <- lapply(result_list, function(x) {
    tmp <- table(unlist(x))
    names(tmp[tmp == 10])
})

selected_features <- list(all_features=all_features, top_features=top_features)
save(selected_features, file='figures/selected_features.RData')


### calculate jaccard index between gene sets ###
# Jaccard index
jaccard_index <- function(x, y){
    intersect <- length(intersect(x, y))
    union <- length(union(x, y))
    return(intersect / union)
}

jaccard_dissimilarity <- function(x, y){
    intersect <- length(intersect(x, y))
    union <- length(union(x, y))
    return(1-(intersect / union))
}

ifelse(dir.exists('figures/jaccard_heatmaps') == F, dir.create('figures/jaccard_heatmaps'))

celltypes <- unique(gsub('.HVG.autosome|.HVG|.autosome|.SLE|.chrX', '', names(selected_features$all_features)))
for(celltype in celltypes){
    tmp <- selected_features$all_features[grep(paste0('^',celltype), names(selected_features$all_features))]
    combinations <- combn(names(tmp), 2)
    mtx <- matrix(0, ncol=5, nrow=5, dimnames = list(names(tmp), names(tmp)))
    for(i in 1:ncol(combinations)){
        mtx[combinations[1,i], combinations[2,i]] <- jaccard_index(tmp[[combinations[1,i]]], tmp[[combinations[2,i]]])
        mtx[combinations[2,i], combinations[1,i]] <- mtx[combinations[1,i], combinations[2,i]]
    }
    diag(mtx) <- 1
    colnames(mtx) <- gsub(celltype, '', colnames(mtx)) %>% gsub('^\\.', '', .)
    rownames(mtx) <- gsub(celltype, '', rownames(mtx)) %>% gsub('^\\.', '', .)
    pdf(paste0('figures/jaccard_heatmaps/', celltype, '.pdf'))
    col = circlize::colorRamp2(c(0, 1), c("white", "red"))
    print(Heatmap(mtx, col=col, name = 'Jaccard Index', column_title=replace.names('CD16+.NK.cells'),
    cluster_columns=FALSE, cluster_rows=FALSE))
    dev.off()
}

result_list <- list()
for(celltype in celltypes){
    result <- lapply(c('.HVG', '.SLE'), function(x){
    jaccard_index(selected_features$all_features[[paste0(celltype, x)]], selected_features$all_features[[paste0(celltype, '.chrX')]])
    })
    names(result) <- c('HVG', 'SLE')
    tmp <- do.call(rbind, result)
    colnames(tmp) <- celltype
    result_list[[celltype]] <- tmp
}
jaccard_matrix <- do.call(cbind, result_list)
colnames(jaccard_matrix) <- replace.names(colnames(jaccard_matrix))

pdf('figures/chrX.geneset.jaccard.heatmap.pdf')
col = circlize::colorRamp2(c(0, 1), c("white", "red"))
Heatmap(t(jaccard_matrix), col=col, name = 'Jaccard Index',
cluster_columns=FALSE, cluster_rows=FALSE)
dev.off()

intersect_list <- list()
for(celltype in celltypes){
    result <- unlist(lapply(c('.HVG', '.SLE'), function(x){
    intersect(selected_features$all_features[[paste0(celltype, x)]], selected_features$all_features[[paste0(celltype, '.chrX')]])
    }))
    intersect_list[[celltype]] <- result
}
sort(table(unlist(intersect_list)))

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
edgeR <- deg.list('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/lupus_Chun/differential.expression/edgeR', 
filter=F)
names(edgeR) <- gsub('_', '.', names(edgeR))

degs <- lapply(celltypes, function(x){
    subset(edgeR[[x]], gene %in% intersect_list[[x]])
})
names(degs) <- celltypes

combined_degs <- bind_rows(degs, .id='celltype')
degs_mtx <- reshape2::dcast(combined_degs, gene ~ celltype, value.var='logFC')
rownames(degs_mtx) <- degs_mtx$gene
degs_mtx <- degs_mtx[,-1]
degs_mtx[is.na(degs_mtx)] <- 0
colnames(degs_mtx) <- replace.names(colnames(degs_mtx))

colnames(degs_mtx) <- replace.names(colnames(degs_mtx))

fdr_matrix <- reshape2::dcast(combined_degs, gene ~ celltype, value.var='FDR')
rownames(fdr_matrix) <- fdr_matrix$gene
colnames(fdr_matrix) <- replace.names(colnames(fdr_matrix))
fdr_matrix <- fdr_matrix[,-1]
significance_matrix <- ifelse(fdr_matrix < 0.05, "*", "ns")
significance_matrix <- ifelse(fdr_matrix < 0.01, "**", significance_matrix)
significance_matrix <- ifelse(fdr_matrix < 0.001, "***", significance_matrix)
significance_matrix[is.na(significance_matrix)] <- ""

pdf('figures/shared_features_deg_heatmap.pdf')
col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
Heatmap(as.matrix(degs_mtx), col=col, name = 'logFC',
        cluster_columns=TRUE, cluster_rows=TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(significance_matrix[i, j], x, y, gp = gpar(fontsize = 10))
        })
dev.off()

### Calculate Jaccard index between chrX features ###
chrX_features <- selected_features$all_features[grep('.chrX', names(selected_features$all_features))]
names(chrX_features) <- gsub('.chrX', '', names(chrX_features))

top_celltypes <- c('CD16+.NK.cells', 'Memory.B.cells','Tcm.Naive.helper.T.cells', 'Regulatory.T.cells')
chrX_features <- chrX_features[top_celltypes]

combinations <- combn(names(chrX_features), 2)
jaccard_mtx <- matrix(0, ncol=4, nrow=4, dimnames = list(names(chrX_features), names(chrX_features)))
for(i in 1:ncol(combinations)){
    jaccard_mtx[combinations[1, i], combinations[2, i]] <- jaccard_index(chrX_features[[combinations[1, i]]], chrX_features[[combinations[2, i]]])
    jaccard_mtx[combinations[2, i], combinations[1, i]] <- jaccard_mtx[combinations[1, i], combinations[2, i]]
}
diag(jaccard_mtx) <- 1

colnames(jaccard_mtx) <- replace.names(colnames(jaccard_mtx))
rownames(jaccard_mtx) <- replace.names(rownames(jaccard_mtx))
jaccard_mtx[lower.tri(jaccard_mtx)] <- 0

pdf('figures/chrX_jaccard_heatmap.pdf')
col = circlize::colorRamp2(c(0, 1), c("white", "red"))
Heatmap(jaccard_mtx, col=col, name = 'Jaccard Index')
dev.off()

edgeR <- deg.list('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/lupus_Chun/differential.expression/edgeR', logfc=0.1)
names(edgeR) <- gsub('_', '.', names(edgeR))
degs <- lapply(top_celltypes, function(x){
    subset(edgeR[[x]], gene %in% chrX_features[[x]])
})
names(degs) <- top_celltypes

combined_degs <- bind_rows(degs, .id='celltype')
degs_mtx <- reshape2::dcast(combined_degs, gene ~ celltype, value.var='logFC')
rownames(degs_mtx) <- degs_mtx$gene
degs_mtx <- degs_mtx[,-1]
degs_mtx[is.na(degs_mtx)] <- 0
colnames(degs_mtx) <- replace.names(colnames(degs_mtx))

fdr_matrix <- reshape2::dcast(combined_degs, gene ~ celltype, value.var='FDR')
rownames(fdr_matrix) <- fdr_matrix$gene
colnames(fdr_matrix) <- replace.names(colnames(fdr_matrix))
fdr_matrix <- fdr_matrix[,-1]
significance_matrix <- ifelse(fdr_matrix < 0.05, "*", "ns")
significance_matrix <- ifelse(fdr_matrix < 0.01, "**", significance_matrix)
significance_matrix <- ifelse(fdr_matrix < 0.001, "***", significance_matrix)
significance_matrix[is.na(significance_matrix)] <- ""

X.immune <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/X.immune.txt', header=FALSE)$V1
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')
katsir <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/Katsir.escape.txt')$Gene.Symbol
SLE <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/SLE.tsv')$Gene
chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name

pdf('figures/chrX_deg_heatmap.pdf')
col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
ann <- rowAnnotation(foo = anno_mark(at = which(rownames(degs_mtx) %in% gene_label), 
                                    labels = rownames(degs_mtx)[rownames(degs_mtx) %in% gene_label]))
Heatmap(as.matrix(degs_mtx), col=col, name = 'logFC',
        cluster_columns=TRUE, cluster_rows=TRUE,
        show_row_names=TRUE, row_names_gp = gpar(fontsize = 5),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(significance_matrix[i, j], x, y, gp = gpar(fontsize = 10))
        })
dev.off()

foo <- subset(average_ensemble_metrics, gene.set=='chrX' & round(MCC, 0) >= 0.7)[,c('celltype','MCC', 'MCC_lower', 'MCC_upper')]
foo$MCC <- round(foo$MCC, 2)
foo$MCC_lower <- round(foo$MCC_lower, 2)
foo$MCC_upper <- round(foo$MCC_upper, 2)
foo

subset(combined_degs, gene %in% c('RPS4X', 'RPL39', 'RPL36A'))


### Considently selected features ###
chrX_features <- selected_features$top_features[grep('.chrX', names(selected_features$top_features))]
names(chrX_features) <- gsub('.chrX', '', names(chrX_features))

top_celltypes <- c('CD16+.NK.cells', 'Memory.B.cells','Tcm.Naive.helper.T.cells', 'Regulatory.T.cells')
chrX_features <- chrX_features[top_celltypes]

edgeR <- deg.list('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/lupus_Chun/differential.expression/edgeR', filter=F)
names(edgeR) <- gsub('_', '.', names(edgeR))
degs <- lapply(top_celltypes, function(x){
    subset(edgeR[[x]], gene %in% chrX_features[[x]])[,c('gene', 'logFC', 'FDR')]
})
names(degs) <- replace.names(top_celltypes)

combined_degs <- bind_rows(degs, .id='celltype')
write.csv(combined_degs, 'figures/top_chrX.consistent.csv', row.names=FALSE)

# Upset plot
chrX_features_mtx <- fromList(chrX_features)
rownames(chrX_features_mtx) <- unique(unlist(chrX_features))
sort(rowSums(chrX_features_mtx), decreasing=FALSE)

pdf('figures/chrX_features_UpSet.pdf', onefile=F)
upset(chrX_features_mtx, order.by = "freq", main.bar.color = "black", sets.bar.color = "black", matrix.color = "black", nset=4)
dev.off()

library(enrichR)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("Enrichr") # Human genes   
}
if (websiteLive) dbs <- listEnrichrDbs()
if (websiteLive) head(dbs)
dbs <- c("GO_Biological_Process_2023", "Reactome_2022", "MSigDB_Hallmark_2020")
if (websiteLive) {
    enriched <- enrichr(chrX_features[[2]], dbs)
}

lapply(enriched, function(x) subset(x, Adjusted.P.value < 0.05))











degs_genes_mts <- degs_genes_mtx[order(rowSums(degs_genes_mtx), decreasing=TRUE),]

rownames(chrX_features_mtx)[rownames(chrX_features_mtx) %in% rownames(escape)]



feature_names <- names(chrX_features)
n <- length(feature_names)
# Initialize a square matrix to store the Jaccard index scores
jaccard_matrix <- matrix(0, n, n, dimnames = list(feature_names, feature_names))

combinations <- combn(feature_names, 2)
mtx <- matrix(0, ncol=22, nrow=22, dimnames = list(feature_names, feature_names))
for(i in 1:ncol(combinations)){
    mtx[combinations[1,i], combinations[2,i]] <- jaccard_index(chrX_features[[combinations[1,i]]], chrX_features[[combinations[2,i]]])
    mtx[combinations[2,i], combinations[1,i]] <- mtx[combinations[1,i], combinations[2,i]]
}
diag(mtx) <- 1

colnames(mtx) <- replace.names(gsub('.chrX', '', colnames(mtx)))
rownames(mtx) <- replace.names(gsub('.chrX', '', rownames(mtx)))

pdf('figures/chrX_jaccard_heatmap.pdf')
col = circlize::colorRamp2(c(0, 1), c("white", "red"))
Heatmap(mtx, col=col, name = 'Jaccard Index')
dev.off()

library(UpSetR)
mtx <- fromList(chrX_features)
colnames(mtx) <- replace.names(gsub('.chrX', '', colnames(mtx)))
rownames(mtx) <- unique(unlist(chrX_features))
pdf('figures/upset_plot.pdf')
upset(mtx, order.by = "freq", main.bar.color = "black", sets.bar.color = "black", matrix.color = "black", nset=16)
dev.off()

X.immune <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/X.immune.txt', header=FALSE)$V1
genes <- rownames(mtx)[rownames(mtx) %in% rownames(escape)]

genes <- names(sort(rowSums(mtx[genes,]), decreasing=TRUE)[1:10])


pdf('figures/chrX_features_heatmap.pdf')
ha = rowAnnotation(foo = anno_mark(at = which(rownames(mtx) %in% genes), 
    labels = genes))
col = circlize::colorRamp2(c(0, 1), c("white", "red"))
Heatmap(as.matrix(mtx), col=col, right_annotation = ha, 
show_row_names = FALSE, show_heatmap_legend = FALSE)
dev.off()

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


load('figures/selected_features.RData')
load('/directflow/SCCGGroupShare/projects/lacgra/eqtls_results.Rdata')
eqtl_features <- lapply(selected_features[['all_features']], function(x){
    x[x %in% sbeqtls$geneid]
})
eqtl_features <- eqtl_features[!sapply(eqtl_features, function(x) length(x) == 0)]

names(unlist(lapply(eqtl_features, function(x) x[x %in% 'IL2RA'])))
subset(female, female$geneid == 'IL2RA')

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
degs <- deg.list('../differential.expression/edgeR', logfc=0.1)
lapply(degs, function(x) subset(x, gene %in% eqtls)[,c('gene', 'logFC', 'FDR')])


eqtls <- unique(unlist(eqtl_features))

subset(female, geneid %in% eqtls)

gwas <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/SLE_GWAS.tsv')

unique(gwas$Gene)[unique(gwas$Gene) %in% eqtls]