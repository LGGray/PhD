library(dplyr)

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')

# Read in ML metrics and combine into a single data frame
# Read in metrics files containing chrX and metrics in filename
metric.files <- list.files('psuedobulk/metrics/', pattern=c('_metrics.*chrX.csv'), full.names=T)
metrics.list <- lapply(metric.files, read.csv)

metric.files <- gsub('_metrics|.csv', '', basename(metric.files))
names(metrics.list) <- metric.files
metrics <- bind_rows(metrics.list, .id='model')

# Write data frame file
write.table(metrics, 'psuedobulk/metrics/chrX.metrics.combined.txt', row.names=FALSE, quote=F, sep='\t')

# Read in ML metrics and combine into a single data frame
# Read in metrics files containing HVG and metrics in filename
metric.files <- list.files('psuedobulk/metrics/', pattern=c('_metrics.*HVG.csv'), full.names=T)
metrics.list <- lapply(metric.files, read.csv)

metric.files <- gsub('_metrics|.csv', '', basename(metric.files))
names(metrics.list) <- metric.files
metrics <- bind_rows(metrics.list, .id='model')

# Write data frame file
write.table(metrics, 'psuedobulk/metrics/HVG.metrics.combined.txt', row.names=FALSE, quote=F, sep='\t')

# # Read in feature counts and add length to metrics data frame
# feature.files <- list.files('ML.models/features/', pattern=c('.txt'), full.names=TRUE)
# feature.list <- lapply(feature.files, function(x) read.csv(x)$Features)
# feature.files <- gsub('_model|.txt', '', basename(feature.files))
# names(feature.list) <- feature.files
# feature.files <- feature.files[feature.files %in% metric.files]
# feature.files <- feature.files[match(feature.files, metric.files)]
# feature.list <- feature.list[feature.files]

# metrics$nFeatures <- unlist(lapply(feature.list, length))

# # Count how many features are in chrX
# metrics$nchrX <- unlist(lapply(feature.list, function(x){
#     sum(x$Features %in% rownames(chrX))
#     }
# ))

# # Write data frame file
# write.table(metrics, 'exp.matrix/metrics/Metrics.combined.txt', row.names=FALSE, quote=F, sep='\t')

# # Read in confusion matrices and save as Rdata
# confusion.files <- list.files('exp.matrix/metrics', pattern=c('confusion'), full.names=T)
# confusion.list <- lapply(confusion.files, function(x) read.csv(x))
# confusion.files <- gsub('_confusion|.csv', '', basename(confusion.files))
# names(confusion.list) <- confusion.files
# save(confusion.list, file='exp.matrix/metrics/confusion.matrix.Rdata')

# metrics.flt[,c(1,3,4,5)]

# features <- lapply(split(metrics.flt, metrics.flt$celltype), function(x){
#     lst <- feature.list[x$model]
#     names(table(unlist(lst)))
# })
