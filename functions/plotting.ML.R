library(ggplot2)
library(dplyr)
library(UpSetR)

load('../datasets/XCI/escapees.Rdata')

setwd('~/external/ClusterScratch/autoimmune.datasets/')
metrics <- read.delim('UC_GSE182270/exp.matrix/metrics/Metrics.combined.txt')

metrics <- metrics %>%
  mutate(clf=gsub('_.+', '', model),
         celltype=gsub('.+_|.chrX', '', model))
pdf('UC_GSE182270/model.comparison.pdf', width = 15)
ggplot(metrics, aes(x=celltype, y=F1, fill=clf)) +
  geom_col(position = 'dodge') +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_fill_discrete(name='Model') +
  xlab('Cell type')
dev.off()

ggplot(metrics, aes(x=celltype, y=F1, colour=clf, size=nFeatures)) +
  geom_point()

metrics.split <- split(metrics, metrics$celltype)
lapply(metrics.split, function(x) x[order(x$F1, decreasing = T),][1,]) %>%
  bind_rows() %>%
  data.frame() -> best.model
best.model

files = paste0('UC_GSE182270/ML.models/features/', best.model$clf, '_model_', best.model$celltype, '.chrX.txt')

#files <- list.files('UC_GSE125527/ML.models/features/', full.names = T)
feature.list <- lapply(files, function(x) read.delim(x)$Features)
names(feature.list) <- gsub('chrX.txt' , '', basename(files))
pdf('UC_GSE182270/best.chrX.models.pdf')
upset(fromList(feature.list), nsets = length(feature.list))
dev.off()

df <- fromList(feature.list)
rownames(df) <- unique(unlist(feature.list))
df$sum <- rowSums(df)
View(df)
hits <- names(which(sort(table(unlist(feature.list))) > 16))

hits %in% rownames(escape)