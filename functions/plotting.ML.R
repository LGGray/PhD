library(ggplot2)
library(dplyr)
library(UpSetR)

# Read in metrics file
metrics <- read.delim('exp.matrix/metrics/Metrics.combined.txt')

# Create new columns for plotting
metrics <- metrics %>%
  mutate(clf=gsub('_.+', '', model),
         celltype=gsub('.+_|.chrX', '', model))
# Plot column graph of F1 score for each model
pdf('model.comparison.pdf', width = 15)
ggplot(metrics, aes(x=celltype, y=F1, fill=clf)) +
  geom_col(position = 'dodge') +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_fill_discrete(name='Model') +
  xlab('Cell type')
dev.off()
# Plot boxplot of nFeatures for each model
pdf('ML.model.nfeatures.pdf', width=15)
ggplot(metrics, aes(x=celltype, y=nFeatures)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=clf)) +
  scale_color_discrete(name='Model') +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  xlab('Cell type')
dev.off()

# Split by cell type and select the best model for each cell type
metrics.split <- split(metrics, metrics$celltype)
lapply(metrics.split, function(x) x[order(x$F1, decreasing = T),][1,]) %>%
  bind_rows() %>%
  data.frame() -> best.model
best.model
# Read in the features for the best models and create upset plot
files = paste0('ML.models/features/', best.model$clf, '_model_', best.model$celltype, '.chrX.txt')
feature.list <- lapply(files, function(x) read.delim(x)$Features)
names(feature.list) <- gsub('chrX.txt' , '', basename(files))
pdf('best.chrX.models.pdf')
upset(fromList(feature.list), nsets = length(feature.list))
dev.off()

# df <- fromList(feature.list)
# rownames(df) <- unique(unlist(feature.list))
# df$sum <- rowSums(df)
# View(df)
# hits <- names(which(sort(table(unlist(feature.list))) > 12))

# hits %in% rownames(escape)

# edgeR <- edgeR.list('psuedobulk', logfc=0.1)
# for(cell in names(feature.list)){
#   edgeR[[paste0(gsub('.+_model_', '', cell), 'edgeR-LRT')]]
# }
