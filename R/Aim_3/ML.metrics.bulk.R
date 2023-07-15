library(dplyr)
library(ggplot2)

metric.files <- list.files('bulk_RNA/', pattern='metrics', full.names=T)
metrics.list <- lapply(metric.files, read.csv)

metric.files <- gsub('metrics_|.csv', '', basename(metric.files))
names(metrics.list) <- metric.files
metrics <- bind_rows(metrics.list, .id='model')

metrics %>%
    mutate(sample=ifelse(grepl('T_cells', model), 'T_cells', 'B_cells'))

# Write data frame file
write.table(metrics, 'bulk_RNA/Metrics.combined.txt', row.names=FALSE, quote=F, sep='\t')

pdf('ML.plots/F1.forest.bulk.female.pdf')
ggplot(metrics.flt, aes(x=F1, y=celltype, color = features)) +
    geom_point(size = 1, position=position_jitter(width = 0.2, height = 0.2, seed = 42)) +
    geom_errorbarh(
        aes(xmin = F1_lower, xmax = F1_upper),
        height = 0.2,
        position=position_jitter(width = 0.2, height = 0.2, seed = 42)) + 
    theme_bw() +
    xlab("F1 score") + ylab("Cell type") + ggtitle("F1 score for each cell type and model type") +
    labs(color='Model')
dev.off()

