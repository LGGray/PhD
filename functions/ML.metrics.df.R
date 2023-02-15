library(dplyr)

load('../../datasets/XCI/chrX.Rdata')

# Read in ML metrics and combine into a single data frame
metric.files <- list.files('exp.matrix/metrics', pattern=c('metrics'), full.names=T)
metrics.list <- lapply(metric.files, read.csv)

metric.files <- gsub('_metrics|.csv', '', basename(metric.files))
names(metrics.list) <- metric.files
metrics <- bind_rows(metrics.list, .id='model')

# Read in feature counts and add length to metrics data frame
feature.files <- list.files('ML.models/features/', full.names=TRUE)
feature.list <- lapply(feature.files, read.csv)
feature.files <- gsub('_model|.txt', '', basename(feature.files))
feature.files <- feature.files[match(feature.files, metric.files)]

metrics$nFeatures <- unlist(lapply(feature.list, nrow))

# Count how many features are in chrX
metrics$nchrX <- unlist(lapply(feature.list, function(x){
    sum(x$Features %in% rownames(chrX))
    }
))

# Write data frame file
write.table(metrics, 'exp.matrix/metrics/Metrics.combined.txt', row.names=FALSE, quote=F)

# # Read in confusion matrices and save as Rdata
# confusion.files <- list.files('exp.matrix/metrics', pattern=c('confusion'), full.names=T)
# confusion.list <- lapply(confusion.files, function(x) read.csv(x))
# confusion.files <- gsub('_confusion|.csv', '', basename(confusion.files))
# names(confusion.list) <- confusion.files
# save(confusion.list, file='exp.matrix/metrics/confusion.matrix.Rdata')
