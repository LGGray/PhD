library(UpSetR)


files <- list.files('ML.models/features', full.names = TRUE)

# Read in files
feature.list <- lapply(files, function(x){
    read.delim(x)[,1]}
)
names(feature.list) <- gsub('.txt', '', basename(files))

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