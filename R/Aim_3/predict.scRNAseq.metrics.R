metric.files <- list.files(path = "psuedobulk/scRNA", pattern = "csv", full.names = TRUE)

metrics <- lapply(metric.files, function(x) {
  read.csv(x)
})
names(metrics) <- gsub('metrics_|.csv', '', basename(metric.files))

library(stringr)
metrics <- dplyr::bind_rows(metrics, .id='celltype')
metrics$celltype <- gsub('.chrX|.HVG|_.+', '', metrics$celltype)
metrics$features <- stringr::str_extract(basename(metric.files), "chrX|HVG")
metrics$features <- factor(metrics$features, levels=c('chrX', 'HVG'))
metrics$group <- stringr::str_extract(basename(metric.files), "adult|child|all")
metrics$group <- factor(metrics$group, levels=c('adult', 'child', 'all'))

metrics <- subset(metrics, celltype != "HSC.MPP")

write.table(metrics, 'ML.plots/chrX_HVG.metrics.txt', sep='\t', row.names=FALSE)

library(ggplot2)
pdf('ML.plots/chrX_HVG_F1.barplot.pdf')
ggplot(metrics, aes(x=celltype, y=F1, fill=group)) +
    geom_bar(stat='identity', position='dodge') +
    facet_grid(. ~ features) +  # Add facet_grid to show scores for each group
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_flip()
dev.off()


pdf('ML.plots/chrX_HVG_AUPRC.barplot.pdf')
ggplot(metrics, aes(x=celltype, y=AUPRC, fill=group)) +
    geom_bar(stat='identity', position='dodge') +
    facet_grid(. ~ features) +  # Add facet_grid to show scores for each group
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_flip()
dev.off()

pdf('ML.plots/chrX_HVG_Kappa.barplot.pdf')
ggplot(metrics, aes(x=celltype, y=Kappa, fill=group)) +
    geom_bar(stat='identity', position='dodge') +
    facet_grid(. ~ features) +  # Add facet_grid to show scores for each group
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_flip()
dev.off()
