library(dplyr)
library(tidyr)
library(stringr)
library(igraph)
library(ggplot2)
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')

gene.set.colours <- c('chrX'='#8A0798', 'autosome'='#7DB176', 'HVG'='#44ABD9', 'HVG.autosome'='#1007D9', 'SLE'='#D90750')

feature.files <- unlist(lapply(1:10, function(split){
  list.files(paste0('split_', split, '/features'), pattern='intersection', full.names=TRUE)
}))
features <- lapply(feature.files, function(x) read.csv(x)$Feature)
names(features) <- gsub('/features/intersection|.csv', '', feature.files)

vif_list <- list()
for(file in unique(gsub('split_.+_', '', names(features)))){
  tbx <- table(unlist(features[grep(paste0('_',file), names(features))]))
  genes <- names(tbx[tbx == 10])
  if(length(genes) < 3){
    vif_list[[file]] <- NA
    next
  } 
  mtx <- readRDS(paste0(file, '.RDS'))
  mtx$ancestry <- ifelse(mtx$ancestry == 'Asian', 1, 0)
  mtx <- mtx[, genes]
  vif_results <- data.frame()
  for(i in colnames(mtx)) {
      formula <- as.formula(paste(i, "~ ."))
      vif_value <- vif(glm.nb(formula, data = mtx))
      vif_results <- rbind(vif_results, data.frame(Feature = i, VIF = mean(vif_value)))
  }
  vif_list[[file]] <- vif_results
}

vif_list <- vif_list[!is.na(vif_list)]
vif_df <- bind_rows(vif_list, .id='celltype') %>%
    mutate(
        gene.set=str_extract(celltype, "HVG.autosome|chrX|autosome|HVG|SLE"),
        celltype=gsub('.HVG.autosome', '', celltype),
        celltype=gsub('.SLE|.chrX|.autosome|.HVG', '', celltype)
      )
vif_df$gene.set <- factor(avg_vif$gene.set, levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))
vif_df$celltype <- replace.names(vif_df$celltype)

pdf('figures/vif_by_celltype.pdf', width=11, height=5)
ggplot(vif_df, aes(x=gene.set, y=VIF, colour=gene.set)) +
    geom_violin() +
    theme_minimal() +
    labs(x='Gene Set', y='VIF') +
    theme(axis.text.x=element_blank(),
    strip.text = element_text(size = 8)) +
    scale_colour_manual(values=gene.set.colours) +
    facet_wrap(~celltype, scales='free')
dev.off()

pdf('figures/vif_by_geneset.pdf')
ggplot(vif_df, aes(x=gene.set, y=VIF, colour=gene.set)) +
    geom_violin() +
    theme_minimal() +
    labs(x='Gene Set', y='VIF') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours)
dev.off()



avg_vif <- lapply(vif_list, function(x) data.frame(mean_VIF=mean(x$VIF)))
avg_vif <- bind_rows(avg_vif, .id='celltype') %>%
    mutate(
        gene.set=str_extract(celltype, "HVG.autosome|chrX|autosome|HVG|SLE"),
        celltype=gsub('.HVG.autosome', '', celltype),
        celltype=gsub('.SLE|.chrX|.autosome|.HVG', '', celltype)
      )

avg_vif$gene.set <- factor(avg_vif$gene.set, levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))






## Testing for correlation between selected features
split_list <- list()
for (split in 1:10){
  feature.files <- list.files(paste0('split_', split, '/features'), pattern='intersection', full.names=TRUE)
  features <- lapply(feature.files, function(x) read.csv(x)$Feature)
  names(features) <- gsub('intersection_|.csv', '', basename(feature.files))
  features <- features[sapply(features, function(x) length(x) > 1)]

  correlation_list <- list()
  for(cell in 1:length(features)) {
      mtx <- readRDS(paste0(names(features)[cell], '.RDS'))
      mtx$class <- ifelse(mtx$class == 'SLE', 1, 0)
      mtx$ancestry <- ifelse(mtx$ancestry == 'Asian', 1, 0)
      mtx <- mtx[, features[[cell]]]

      vif_results <- data.frame()
      for (feature in colnames(mtx)) {
          formula <- as.formula(paste(feature, "~ ."))
          vif_value <- vif(glm.nb(formula, data = mtx))
          vif_results <- rbind(vif_results, data.frame(Feature = feature, VIF = mean(vif_value)))
      }


      mtx_scaled <- data.frame(cbind(class=mtx$class, mtx_scaled))
      car::vif(lm(HLA.A ~ HLA.DQA1 + HLA.DQB1, data = mtx_scaled))


      pca <- prcomp(mtx, center=TRUE, scale=TRUE)

      

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

      # Create the graph object
      graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
      assortativity_value <- assortativity_degree(graph, directed = FALSE)

      # print(paste0(names(features)[i], ': ', assortativity_value))

      tmp <- data.frame(celltype=names(features)[cell], assortativity=assortativity_value)
      correlation_list[[names(features)[cell]]] <- tmp
  }
  correlation_list <- dplyr::bind_rows(correlation_list, .id=NULL)
  split_list[[split]] <- correlation_list
}

correlation_df <- dplyr::bind_rows(split_list, .id='split')
correlation_df <- correlation_df %>%
    mutate(
        gene.set=str_extract(celltype, "HVG.autosome|chrX|autosome|HVG|SLE"),
        celltype=gsub('.HVG.autosome', '', celltype),
        celltype=gsub('.SLE|.chrX|.autosome|.HVG', '', celltype) 
    )
correlation_df$gene.set <- factor(correlation_df$gene.set, levels=c('chrX', 'autosome', 'HVG', 'HVG.autosome', 'SLE'))
correlation_df[is.na(correlation_df$assortativity),'assortativity'] <- 0
write.csv(correlation_df, 'figures/assortativity.csv')

pdf('figures/assortativity.pdf')
ggplot(correlation_df, aes(x=gene.set, y=assortativity, colour=gene.set)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='Gene Set', y='Number of Features') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours)
dev.off()

avg_correlation_df <- correlation_df %>%
    group_by(celltype, gene.set) %>%
    summarise(mean=mean(assortativity), sd=sd(assortativity))

pdf('figures/avg_assortativity.pdf')
ggplot(avg_correlation_df, aes(x=gene.set, y=mean, colour=gene.set)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='Gene Set', y='Number of Features') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours)
dev.off()