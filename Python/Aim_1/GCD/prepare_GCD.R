library(Seurat)
library(igraph)
library(edgeR)

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/Seurat2PB.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')

chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name

if (!dir.exists('GCD')) {
  dir.create('GCD')
}

load('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/Aim_1/combined_fdr_list.Rdata')

pbmc <- readRDS(commandArgs(trailingOnly = TRUE)[1])

for(cell in names(combined_fdr_list)) {
    if(! replace.names(gsub('_', '.', cell)) %in% levels(pbmc))
        next
    control <- subset(pbmc, cellTypist == replace.names(gsub('_', '.', cell)) & condition == 'control')
    disease <- subset(pbmc, cellTypist == replace.names(gsub('_', '.', cell)) & condition == 'disease')

    if(length(unique(control$individual)) < 2 | length(unique(disease$individual)) < 2) {
        next
    }

    if(!dir.exists(paste('GCD/', cell))) {dir.create(paste0('GCD/', cell))}
    # Create pseudobulked expression matrix - Control
    control_PB <- Seurat2PB(control, sample='individual', cluster='individual', assay='RNA')$counts
    # Remove genes with zero counts
    control_PB <- control_PB[rowSums(control_PB) > 0,]
    # Filter for X chrommosome genes
    control_PB <- control_PB[rownames(control_PB) %in% chrX,]

    # Create pseudobulked expression matrix - Disease
    disease_PB <- Seurat2PB(disease, sample='individual', cluster='individual', assay='RNA')$counts
    # Remove genes with zero counts
    disease_PB <- disease_PB[rowSums(disease_PB) > 0,]
    # Filter for X chrommosome genes
    disease_PB <- disease_PB[rownames(disease_PB) %in% chrX,]

    # Find common genes between control and disease
    common_genes <- intersect(rownames(control_PB), rownames(disease_PB))

    # Filter for common genes
    control_PB <- control_PB[common_genes,]
    disease_PB <- disease_PB[common_genes,]

    # Calculate spearman correlation
    control_coexpression <- cor(t(control_PB), method='spearman')
    # Create igraph object
    control_g <- graph_from_adjacency_matrix(control_coexpression, mode='undirected', weighted=TRUE, diag=FALSE)
    # Save graph
    write.graph(control_g, file=paste0('GCD/', cell, '/control.gw'), format="leda")


    # Calculate spearman correlation
    disease_coexpression <- cor(t(disease_PB), method='spearman')
    # Create igraph object
    disease_g <- graph_from_adjacency_matrix(disease_coexpression, mode='undirected', weighted=TRUE, diag=FALSE)
    # Save graph
    write.graph(disease_g, file=paste0('GCD/', cell, '/disease.gw'), format="leda")
}


# # plot both graphs 
# library(ComplexHeatmap)
# library(circlize)

# # Load the graph
# ht1 <- Heatmap(as.matrix(control_coexpression), name = "Rho",
# column_title = "Control", 
# show_row_names = FALSE, show_column_names = FALSE, 
# cluster_rows = TRUE, cluster_columns = TRUE, 
# col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
# ht2 <- Heatmap(as.matrix(disease_coexpression), name = "Rho",
# column_title = "Disease",
# show_row_names = FALSE, show_column_names = FALSE,
# cluster_rows = TRUE, cluster_columns = TRUE,
# col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

# pdf('GCD/control_heatmap.pdf')
# ht1
# dev.off()

# pdf('GCD/disease_heatmap.pdf')
# ht2
# dev.off()

# # Compute the difference matrix
# difference_matrix <- disease_coexpression - control_coexpression

# # Sum the absolute differences for each node
# total_differences <- apply(abs(difference_matrix), 1, sum)

# # Rank the nodes based on the total differences
# node_ranks <- sort(total_differences)