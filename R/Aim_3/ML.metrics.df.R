library(dplyr)

load('../../datasets/XCI/chrX.Rdata')

# Read in ML metrics and combine into a single data frame
# Read in metrics files containing chrX and metrics in filename
metric.files <- list.files('exp.matrix/metrics/', pattern=c('metrics', 'chrX', 'HVG'), full.names=T)
metrics.list <- lapply(metric.files, read.csv)

metric.files <- gsub('_metrics|.csv', '', basename(metric.files))
names(metrics.list) <- metric.files
metrics <- bind_rows(metrics.list, .id='model')

# Read in feature counts and add length to metrics data frame
feature.files <- list.files('ML.models/features/', pattern=c('.txt'), full.names=TRUE)
feature.list <- lapply(feature.files, read.csv)
feature.files <- gsub('_model|.txt', '', basename(feature.files))
names(feature.list) <- feature.files
feature.files <- feature.files[match(feature.files, metric.files)]
feature.list <- feature.list[feature.files]

metrics$nFeatures <- unlist(lapply(feature.list, nrow))

# # Count how many features are in chrX
# metrics$nchrX <- unlist(lapply(feature.list, function(x){
#     sum(x$Features %in% rownames(chrX))
#     }
# ))

# Write data frame file
write.table(metrics, 'exp.matrix/metrics/Metrics.combined.txt', row.names=FALSE, quote=F, sep='\t')

# # Read in confusion matrices and save as Rdata
# confusion.files <- list.files('exp.matrix/metrics', pattern=c('confusion'), full.names=T)
# confusion.list <- lapply(confusion.files, function(x) read.csv(x))
# confusion.files <- gsub('_confusion|.csv', '', basename(confusion.files))
# names(confusion.list) <- confusion.files
# save(confusion.list, file='exp.matrix/metrics/confusion.matrix.Rdata')


# metrics.flt <- metrics %>% 
#     filter(F1 >= 0.8) %>%
#     mutate(celltype = gsub('.+_', '', model)) %>%
#     arrange(celltype)

# metrics.flt[,c(1,3,4,5)]

# features <- lapply(split(metrics.flt, metrics.flt$celltype), function(x){
#     lst <- feature.list[x$model]
#     names(table(unlist(lst)))
# })
