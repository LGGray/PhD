library(dplyr)

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')

# Read in ML metrics and combine into a single data frame
metric.files <- list.files('psuedobulk/metrics/', pattern=c('_metrics.*.csv'), full.names=T)
metrics.list <- lapply(metric.files, read.csv)

metric.files <- gsub('_metrics|.csv', '', basename(metric.files))
names(metrics.list) <- metric.files
metrics <- bind_rows(metrics.list, .id='model')

# Write data frame file
write.table(metrics, 'psuedobulk/metrics/metrics.combined.txt', row.names=FALSE, quote=F, sep='\t')
