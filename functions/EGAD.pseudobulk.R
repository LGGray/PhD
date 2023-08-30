library(EGAD)
library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)

load('/directflow/SCCGGroupShare/projects/lacgra/CoExpNets/bin/run_GBA.Rdata')
load('../../CoExpNets/bin/run_GBA.Rdata')
load('../../CoExpNets/bin/helper_functions.r')

# Load data and psuedobulk
pbmc <- readRDS('pbmc.female.control-managed.RDS')

cell = 'Non-classical monocytes'
pbmc.cell <- subset(pbmc, cellTypist == cell)

# Keep genes with expression in 5% of cells
keep <- rowSums(pbmc.cell@assays$RNA@counts > 0) > ncol(pbmc.cell) * 0.05
features <- names(keep[keep == T])
pbmc.cell <- subset(pbmc.cell, features=features)

# Psudobulking by summing counts
expr <- AggregateExpression(pbmc.cell, group.by='individual', slot='counts')$RNA
expr <- expr[,(colSums(expr) > 0)]

# Calculate network edges using correlations and rank standardize
network = EGAD::build_coexp_network(expr)

# plot heatmap of correlation matrix
pdf('EGAD.heatmap.cor.pdf')
heatmap.2(network, dendrogram='none', trace='none', col=rev(colorRampPalette(c('blue', 'white', 'red'))(100)), scale='none')
dev.off()