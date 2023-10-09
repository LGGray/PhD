library(ggplot2)
library(ggrepel)
library(dplyr)

metrics.path <- list.files('psuedobulk/ML.models/ensemble', pattern='metrics_', full.names=TRUE)

metrics <- lapply(metrics.path, function(x) {
  df <- read.csv(x)
  df <- cbind(celltype=gsub('.csv', '', basename(x)), df)
  df
})
metrics <- metrics %>% bind_rows() %>%
    mutate(perm=ifelse(grepl('perm', celltype), 'post', 'pre')) %>%
    mutate(celltype=gsub('perm_metrics_|metrics_', '', celltype)) %>%
    mutate(perm=factor(perm, levels=c('pre', 'post'))) %>%
    mutate(features=ifelse(grepl('chrX', celltype), 'chrX', 'HVG')) %>%
    mutate(features=factor(features, levels=c('chrX', 'HVG'))) %>%
    mutate(celltype=gsub('.chrX|.HVG', '', celltype)) %>%
    arrange(celltype)

# ggplot boxplot of F1 score on y axis and perm on x axis but split by features
pdf('psuedobulk/ML.plots/feature_perm_f1.pdf')
ggplot(metrics, aes(x=perm, y=F1, fill=features)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(0.2), alpha = 0.5) +
  facet_wrap(~features) +
  ylab('F1 score') + xlab('')
dev.off()

pdf('psuedobulk/ML.plots/feature_perm_auc.pdf')
ggplot(metrics, aes(x=perm, y=AUC, fill=features)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(0.2), alpha = 0.5) +
  facet_wrap(~features) +
  ylab('F1 score') + xlab('')
dev.off()

features.list <- list.files('psuedobulk/ML.models/ensemble/features', pattern='.txt', full.names=TRUE)
features <- lapply(features.list, function(x) {
  tmp <- read.csv(x, sep='\t', header=TRUE)
  data.frame(celltype=gsub('perm.|.chrX.txt|.HVG.txt', '', basename(x)), 
  size=nrow(tmp),
  perm=ifelse(grepl('perm.', x), 'post', 'pre'),
  features=ifelse(grepl('chrX', x), 'chrX', 'HVG'))
})
features <- bind_rows(features)

# Combine all data together
metrics$id <- paste(metrics$celltype, metrics$perm, metrics$features, sep='_')
features$id <- paste(features$celltype, features$perm, features$features, sep='_')
results <- merge(metrics, features, by='id') %>%
  select(-c(id, celltype.y, perm.y, features.y))
colnames(results) <- gsub('.x', '', colnames(results))

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')
results$celltype <- replace.names(results$celltype)

write.table(results, 'psuedobulk/ML.models/ensemble/feature_size_metrics.txt', sep='\t', row.names=FALSE, quote=FALSE)

# Dot plot of feature size on x axis, F1 score on y axis, fill = AUC and size = feature size
pdf('psuedobulk/ML.plots/feature_size_f1_auc.pdf')
ggplot(results, aes(y=celltype, x=F1, color=AUC)) +
    geom_point(aes(size=size)) +
    scale_color_gradient2(low = "white", mid = "blue", high = "red", 
                                             midpoint = 0.5, limits = c(0, 1), name = "AUC") +
    ylab('') + xlab('F1 score') +
    geom_vline(xintercept = 0.8, linetype = 'dotted') +
    facet_wrap(~features+perm) +
    theme(panel.spacing=unit(1,"lines")) +
    scale_x_continuous(limits = c(0, 1.0), expand = c(0, 0), breaks = seq(0, 1.0, by = 0.2))
dev.off()
  
# Scatter plot pre and post feature size, F1 score and AUC
results.split <- split(results, results$features)
chrX.results <- split(results.split[['chrX']], results.split[['chrX']]$perm)
chrX.results <- merge(chrX.results[['pre']], chrX.results[['post']], by='celltype', suffixes=c('.pre', '.post'))

pdf('psuedobulk/ML.plots/f1.chrX.scatter.pdf')
ggplot(chrX.results, aes(x = F1.pre, y = F1.post, color = celltype)) +
geom_point() +
geom_text_repel(aes(label = celltype), size = 3, box.padding = 0.5, color = "black") +
theme(legend.position="none") +
ylab('F1 score post') + xlab('F1 score pre')
dev.off()
pdf('psuedobulk/ML.plots/auc.chrX.scatter.pdf')
ggplot(chrX.results, aes(x = AUC.pre, y = AUC.post, color = celltype)) +
geom_point() +
geom_text_repel(aes(label = celltype), size = 3, box.padding = 0.5, color = "black") +
theme(legend.position="none") +
ylab('AUC post') + xlab('AUC pre')
dev.off()

HVG.results <- split(results.split[['HVG']], results.split[['HVG']]$perm)
HVG.results <- merge(HVG.results[['pre']], HVG.results[['post']], by='celltype', suffixes=c('.pre', '.post'))

pdf('psuedobulk/ML.plots/f1.HVG.scatter.pdf')
ggplot(HVG.results, aes(x = F1.pre, y = F1.post, color = celltype)) +
geom_point() +
geom_text_repel(aes(label = celltype), size = 3, box.padding = 0.5, color = "black") +
theme(legend.position="none") +
ylab('F1 score post') + xlab('F1 score pre')
dev.off()
pdf('psuedobulk/ML.plots/auc.HVG.scatter.pdf')
ggplot(HVG.results, aes(x = AUC.pre, y = AUC.post, color = celltype)) +
geom_point() +
geom_text_repel(aes(label = celltype), size = 3, box.padding = 0.5, color = "black") +
theme(legend.position="none") +
ylab('AUC post') + xlab('AUC pre')
dev.off()


# dot plot pf x=f1, y= cell type, fill = AUC and size =feature size

