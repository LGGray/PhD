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


upset(fromList(chrX.features), nsets=length(chrX.features))

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

