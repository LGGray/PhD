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
    #write.graph(control_g, file=paste0('GCD/', cell, '/control.gw'), format="leda")

    # Calculate spearman correlation
    disease_coexpression <- cor(t(disease_PB), method='spearman')
    # Create igraph object
    disease_g <- graph_from_adjacency_matrix(disease_coexpression, mode='undirected', weighted=TRUE, diag=FALSE)
    # Save graph
    #write.graph(disease_g, file=paste0('GCD/', cell, '/disease.gw'), format="leda")

    control_eigenvector <- eigen_centrality(control_g)$vector
    disease_eigenvector <- eigen_centrality(disease_g)$vector
    eigenvector_difference <- abs(disease_eigenvector - control_eigenvector)
    save(eigenvector_difference, file=paste0('GCD/', cell, '/eigenvector_difference.Rdata'))

    control_degree <- igraph::degree(control_g)
    disease_degree <- igraph::degree(disease_g)
    degree_difference <- disease_degree - control_degree
    save(degree_difference, file=paste0('GCD/', cell, '/degree_difference.Rdata'))
}

eigenvector.files <- list.files('GCD', recursive=TRUE, full.names=TRUE, pattern='eigenvector_difference.Rdata')
eigenvector_difference <- lapply(eigenvector.files, function(x) {
    load(x)
    return(data.frame(gene=names(eigenvector_difference), difference=eigenvector_difference))
})
names(eigenvector_difference) <- replace.names(gsub('_', '.', sapply(strsplit(eigenvector.files, '/'), function(x) x[length(x)-1])))
eigenvector_df <- bind_rows(eigenvector_difference, .id='celltype')

pdf('GCD/eigenvector_difference.pdf')
ggplot(eigenvector_df, aes(x=difference, y=celltype)) + 
geom_boxplot() +
theme_minimal() +
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
    
    
strsplit(dirname(x), '/')[[1]][2]

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

