library(UpSetR)
library(dplyr)
library(ggplot2)

metrics <- read.delim('exp.matrix/metrics/Metrics.combined.txt')
sorted.metrics <- metrics %>%
  mutate(obs=gsub('.+_', '', model)) %>%
  group_by(obs) %>%
  group_split() %>%
  lapply(function(x){
    x %>%
      arrange(desc(F1), desc(AUC))
  }) %>%
  bind_rows() %>%
  data.frame() %>%
  subset(., F1 > 0.7)

res <- lapply(split(sorted.metrics, sorted.metrics$obs), head,1) %>%
  bind_rows() %>%
  data.frame()
nrow(res[grep('.chrX', res$model),][,c(1,5:8)])

# Read in model features
feature.files <- list.files('ML.models/features', full.names = TRUE)
# Read in feature.files
feature.list <- lapply(feature.files, function(x){
    read.delim(x)[,1]}
)
names(feature.list) <- gsub('.txt', '', basename(feature.files))

# Get the different celltypes
cells <- unique(gsub('.+_model_|.deg|.chrX', '', names(feature.list)))
# Split names by '_' and get the third element
cells <- unlist(lapply(names(feature.list), function(x) strsplit(x, '_')[[1]][3]))

# find overlap between the different celltypes
overlap <- lapply(cells, function(x) {
  # Get the names of the features for the celltype
  features <- feature.list[grepl(x, names(feature.list))]
  # Get the features that are in features$features
    Reduce(intersect, features)
  # Get the number of features that are in all the models
  length(overlap)
})

# Find common X-linked features between the different celltypes
chrX.features <- lapply(cells, function(x){
    features <- feature.list[grepl(paste0(x, '.chrX'), names(feature.list))]
    Reduce(intersect, features)
})
names(chrX.features) <- cells

# remove empty lists
chrX.features <- chrX.features[sapply(chrX.features, length) > 0]

pdf('psuedobulk/ML.plots/Selected.features.upset.pdf')
upset(fromList(lapply(features, function(x) x$Features)), nsets=length(features))
dev.off()

# Heatmap of features in each celltype
features.mtx <- fromList(lapply(features, function(x) x$Features))
rownames(features.mtx) <- unique(unlist(lapply(features, function(x) x$Features)))

hits <- c('IGBP1', 'TSC22D3', 'CD40LG', 'BTK', 'TLR7', 'XIST', 'IL2RG', 'TMSB4X', 'CYBB', 'KDM6A', 'USP9X', 'MECP2', 'G6PD', 'SH2D1A') 
ha = rowAnnotation(foo = anno_mark(at = which(rownames(features.mtx) %in% hits), 
    labels = rownames(features.mtx)[which(rownames(features.mtx) %in% hits)]))
pdf('psuedobulk/ML.plots/Selected.features.heatmap.pdf')
Heatmap(as.matrix(features.mtx), show_row_names=FALSE, right_annotation = ha, show_heatmap_legend = FALSE)
dev.off()

df <- fromList(chrX.features)
rownames(df) <- unique(unlist(chrX.features))
df$size = rowSums(df)
df <- df[order(df$size, decreasing = TRUE),]

Reduce(intersect, chrX.features)

p
pdf('best.chrX.models.pdf', width=10, height=10)
upset(fromList(lst), nsets=length(lst))
dev.off()

df <- fromList(lst)
rownames(df) <- unique(unlist(lst))
df$size = rowSums(df)
df <- df[order(df$size, decreasing = TRUE),]
rownames(df)[1:10]

load('../../datasets/XCI/escapees.Rdata')

rownames(df[df$size>12,])[rownames(df[df$size>12,]) %in% rownames(escape)]

edgeR <- edgeR.list('psuedobulk/', filter=F)
names(edgeR) <- gsub('.edgeR-LRT', '', names(edgeR))

for(i in 1:3){
  x <- (colnames(df)[which(df[i,] == 1)])
  print(unique(gsub('.+_model_|\\.chrX', '', x)))
  print('\n')
}



####################
metrics <- read.delim('exp.matrix/metrics/Metrics.combined.txt')

# Filter for F1 > 0.8
metrics.flt <- metrics %>% 
    filter(F1 >= 0.8) %>%
    mutate(celltype = gsub('.+_', '', model)) %>%
    mutate(features = gsub('^.*\\.', '', model)) %>%
    arrange(celltype)

feature.files <- list.files('ML.models/features/', full.names = TRUE)
feature.files <- feature.files[gsub('_model|.txt', '', basename(feature.files)) %in% metrics.flt$model]
feature.list <- lapply(feature.files, function(x) read.csv(x)$Features)
names(feature.list) <- gsub('_model|.txt', '', basename(feature.files))

intersection <- lapply(unique(gsub('.chrX|.HVG', '', metrics.flt$celltype)), function(x){
  chrX <- unique(unlist(feature.list[grep(paste0(x, '.chrX'), names(feature.list))]))
  HVG <- unique(unlist(feature.list[grep(paste0(x, '.HVG'), names(feature.list))]))
  intersect(chrX, HVG)
})
names(intersection) <- unique(gsub('.chrX|.HVG', '', metrics.flt$celltype))