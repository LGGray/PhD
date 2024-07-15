library(dplyr)
library(tidyr)
library(stringr)
library(igraph)
library(ggplot2)

gene.set.colours <- c('chrX'='#8A0798', 'autosome'='#7DB176', 'HVG'='#44ABD9', 'HVG.autosome'='#1007D9', 'SLE'='#D90750')


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
        mtx <- mtx[, features[[cell]]]
        # drop ancestry from mtx
        mtx <- mtx[,colnames(mtx) != 'ancestry']
        correlation <- cor(mtx, method='spearman')

        
        n <- ncol(mtx)
        adj_matrix <- matrix(0, n, n)

        # Loop over each pair of variables
        for (i in 1:(n-1)) {
          for (j in (i+1):n) {
            test_result <- cor.test(mtx[, i], mtx[, j], method="spearman")
            # If the correlation is significant, set to 1
            if (test_result$p.value < 0.05) {
              adj_matrix[i, j] <- 1
              adj_matrix[j, i] <- 1
            }
          }
        }

        # Create the graph object
        graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE, weighted=TRUE)
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