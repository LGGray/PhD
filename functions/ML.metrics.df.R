library(dplyr)

# Read in ML metrics and combine into a single data frame
metric.files <- list.files('exp.matrix/metrics', pattern=c('metrics'))
metrics.list <- lapply(metric.files, read.csv)

metric.files <- gsub('_metrics|.csv', '', metric.files)
names(metrics.list) <- metric.files
metrics <- bind_rows(metrics.list, .id='model')

# Read in feature counts and add length to metrics data frame
feature.files <- list.files('ML.models/features/', full.names=TRUE)
feature.list <- lapply(feature.files, read.csv)
feature.files <- gsub('_model|.txt', '', basename(feature.files))
feature.files <- feature.files[match(feature.files, metric.files)]

metrics$nFeatures <- unlist(lapply(feature.list, nrow))

# Write data frame file
write.table(metrics, 'Metrics.combined.txt', row.names=FALSE, quote=F)

# Read in confusion matrices and save as Rdata
confusion.files <- list.files('exp.matrix/metrics', pattern=c('confusion'))
confusion.list <- lapply(confusion.files, read.csv, row.names=c('Control', 'Disease'), col.names=c('Control', 'Disease'))
confusion.files <- gsub('_confusion|.csv', '', confusion.files)
names(confusion.list) <- confusion.files
