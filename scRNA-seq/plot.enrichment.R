setwd('~/external/ClusterHome/')

library(ggplot2)
library(ggrepel)
library(tidyverse)
source('R_code/functions/edgeR.list.R')
source('R_code/functions/fishers.test.degs.R')
load('datasets/XCI/escapees.Rdata')

dataset <- 'GSE193770'

deg.list <- edgeR.list(paste0('datasets/', dataset, '/edgeR-LRT/'), filter = F)
enrichment <- lapply(deg.list, function(x){
  fishers.test.degs(x, genes=rownames(escape), logfc=0.5)
})

df <- lapply(enrichment, function(x) data.frame(p.value=x$p.value, OR=x$estimate))
df <- dplyr::bind_rows(df, .id='column_label')
df$count <- unlist(lapply(deg.list, function(x){
  nrow(subset(x, FDR < 0.5 & logFC > 0.5 & gene %in% rownames(escape)))
}))

pdf(paste0('datasets/', dataset, '/escape.enrichment.pdf'))
ggplot(aes(x = p.value, y = OR, color = p.value < 0.05 & OR > 0, size=count), data = df) +
  geom_point() +
  geom_text_repel(aes(label = column_label), size=5, color = 'black', data = df %>% filter(abs(log(OR)) > log(1))) + 
  theme_minimal() +
  labs(y = 'odds ratio', title = 'XCI escape enrichment', x = 'probability')
dev.off()


