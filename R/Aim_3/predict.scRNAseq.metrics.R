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

source('../../PhD/functions/replace.names.R')
metrics$celltype <- replace.names(metrics$celltype)

write.table(metrics, 'ML.plots/chrX_HVG.metrics.txt', sep='\t', row.names=FALSE)

library(ggplot2)

pdf('ML.plots/chrX_HVG_F1.dotplot.pdf')
ggplot(metrics, aes(x=celltype, y=F1, fill=group)) +
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8)) +
    geom_errorbar(aes(ymin=F1_lower, ymax=F1_upper), width=0.2, position=position_dodge(0.8)) +  # Add error bars
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "black") +  # Add dotted line at 0.8
    facet_grid(. ~ features) +  # Add facet_grid to show scores for each group
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('') +
    coord_flip()
dev.off()

pdf('ML.plots/chrX_HVG_AUPRC.dotplot.pdf')
ggplot(metrics, aes(x=celltype, y=AUPRC, fill=group)) +
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8)) +
    geom_errorbar(aes(ymin=AUPRC_lower, ymax=AUPRC_upper), width=0.2, position=position_dodge(0.8)) +  # Add error bars
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "black") +  # Add dotted line at 0.8
    facet_grid(. ~ features) +  # Add facet_grid to show scores for each group
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('') +
    coord_flip()
dev.off()


