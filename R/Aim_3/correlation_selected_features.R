library(dplyr)
library(tidyr)
library(stringr)
library(igraph)
library(ggplot2)
library(ggpubr)
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')

gene.set.colours <- c('chrX'='#8A0798', 'autosome'='#7DB176', 'HVG'='#44ABD9', 'HVG.autosome'='#1007D9', 'SLE'='#D90750')

load('figures/selected_features.RData')

assortativity_list <- list()
for(file in names(selected_features$all_features)){
  features <- selected_features$all_features[[file]]
  mtx <- readRDS(paste0(file, '.RDS'))
  if('ancestry' %in% features){
    mtx$ancestry <- ifelse(mtx$ancestry=='Asian', 1, 0)
  }
  mtx <- scale(mtx[,features])

  n <- ncol(mtx)
  adj_matrix <- matrix(0, n, n)
  # Loop over each pair of variables
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      test_result <- cor.test(mtx[, i], mtx[, j], method="spearman")
      # If the correlation is significant, set to 1
      if (test_result$p.value < 0.05) {
        adj_matrix[i, j] <- test_result$estimate
        adj_matrix[j, i] <- test_result$estimate
      } else {
          adj_matrix[i, j] <- 0
          adj_matrix[j, i] <- 0
      }
    }
  }
  graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  assortativity_value <- assortativity_degree(graph, directed = FALSE)
  assortativity_list[[file]] <- assortativity_value
}
assortativity_df <- data.frame(celltype=names(assortativity_list), assortativity=unlist(assortativity_list))
assortativity_df$nFeatures <- unlist(lapply(selected_features$all_features, length))
assortativity_df$gene.set <- factor(str_extract(assortativity_df$celltype, "HVG.autosome|chrX|autosome|HVG|SLE"), levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))
assortativity_df$celltype <- replace.names(gsub(".HVG.autosome|.chrX|.autosome|.HVG|.SLE", '', assortativity_df$celltype))

assortativity_metrics <- merge(assortativity_df, average_ensemble_metrics, by=c('celltype', 'gene.set'))
assortativity_metrics$gene.set <- factor(assortativity_metrics$gene.set, levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))
write.csv(assortativity_metrics, 'figures/assortativity_metrics.csv', row.names=FALSE)


lapply(split(assortativity_metrics, assortativity_metrics$gene.set), function(x) cor.test(x$assortativity, x$MCC, method='spearman'))

pairwise.wilcox.test(assortativity_metrics$assortativity, assortativity_metrics$gene.set, p.adjust.method='fdr')
pairwise.wilcox.test(assortativity_metrics$nFeatures, assortativity_metrics$gene.set, p.adjust.method='fdr')

pdf('figures/assortativity_MCC_nFeatures.pdf', width=10, height=5)
ggplot(assortativity_metrics, aes(x=MCC, y=assortativity, colour=gene.set)) +
    geom_point(aes(size=nFeatures)) +
    geom_smooth(method='lm', se=FALSE, color='black') +
    theme_minimal() +
    scale_colour_manual(values=gene.set.colours) +
    facet_wrap(~gene.set, ncol=5, nrow=1) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", color='black', label.y = -0.4)
dev.off()


pdf('figures/assortativity_geneset.pdf')
ggplot(assortativity_df, aes(x=gene.set, y=assortativity, colour=gene.set)) +
    geom_jitter(width=0.2, aes(size=nFeatures)) +
    geom_boxplot(outlier.shape=NA, color='black', fill=NA) +
    theme_minimal() +
    scale_colour_manual(values=gene.set.colours)
dev.off()


pdf('figures/nFeatures_geneset.pdf')
comparisons <- list(c('chrX', 'autosome'), c('chrX', 'HVG'), c('chrX', 'HVG.autosome'), c('chrX', 'SLE'))
ggplot(assortativity_df, aes(x=gene.set, y=nFeatures, colour=gene.set)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape=NA, color='black', fill=NA) +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, 
    test='wilcox.test', color='black') +
    theme_minimal() +
    scale_colour_manual(values=gene.set.colours)
dev.off()

