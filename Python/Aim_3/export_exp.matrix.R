library(Seurat)
library(edgeR)

# This Rscript reads in a given Seurat object and exports the SCTransformed expression matrix for machine Learning analysis
options(warn=-1)

# Load X chromosome genes
chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name

# Load Pseudobulk function
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/Seurat2PB.R')

# Take in command line arguments i.e the Seurat object filename
file <- commandArgs(trailingOnly=TRUE)[1]

pbmc <- readRDS(file)

###### Pseudobulk analysis ######

# Export psuedobulked expression matrix, subsetted by X chromosome genes for each cell type
# Gene must be expressed in at least 5% of individuals
for(cell in levels(pbmc)){
    pbmc.subset <- subset(pbmc, cellTypist == cell)
    bulk <- Seurat2PB(pbmc.subset, sample='individual', cluster='condition', assay='RNA')$counts
    keep <- apply(bulk, 1, function(x) sum(x > 0) > ncol(bulk) * 0.05)
    bulk <- bulk[keep,]
    bulk <- data.frame(t(bulk[rownames(bulk) %in% chrX,]))
    meta <- unique(pbmc.subset@meta.data[,c('condition', 'individual')])
    meta$cellCount <- sapply(meta$individual, function(id) sum(pbmc.subset$individual == id))
    meta <- meta[match(gsub('_.+', '', rownames(bulk)), meta$individual),]
    bulk <- cbind(class=meta$condition, individual=meta$individual, cellCount=meta$cellCount, bulk)
    cell <- gsub('/| |-', '.', cell)
    saveRDS(bulk, paste0('pseudobulk/', cell, '.chrX.RDS'))
}



# 
#     pbmc.subset <- subset(pbmc, cellTypist == cell)
#     # Pseudobulk by average expression
#     bulk <- AverageExpression(pbmc.subset, slot='counts', group.by='individual')$RNA
#     # Filter genes that are expressed in at least 5% of individuals
#     keep <- apply(bulk, 1, function(x) sum(x > 0) > ncol(bulk) * 0.05)
#     bulk <- bulk[keep,]
#     bulk <- data.frame(t(bulk[rownames(bulk) %in% rownames(chrX),]))
#     meta <- unique(pbmc.subset@meta.data[,c('condition', 'individual')])
#     meta$cellCount <- sapply(meta$individual, function(id) sum(pbmc.subset$individual == id))
#     meta <- meta[match(rownames(bulk), meta$individual),]
#     bulk <- cbind(class=meta$condition, individual=meta$individual, cellCount=meta$cellCount, bulk)
#     cell <- gsub('/| |-', '.', cell)
#     saveRDS(bulk, paste0('psuedobulk/', cell, '.chrX.RDS'))
# }

# Export the pseudobulked expression matrix, subsetted by HVG for each cell type
for(cell in levels(pbmc)){
    pbmc.subset <- subset(pbmc, cellTypist == cell)
    pbmc.subset <- FindVariableFeatures(pbmc.subset, nfeatures=2000)
    bulk <- Seurat2PB(pbmc.subset, sample='individual', cluster='condition', assay='RNA')$counts
    # Filter genes that are expressed in at least 5% of individuals
    keep <- apply(bulk, 1, function(x) sum(x > 0) > ncol(bulk) * 0.05)
    bulk <- bulk[keep,]
    bulk <- data.frame(t(bulk[rownames(bulk) %in% VariableFeatures(pbmc.subset),]))
    
    meta <- unique(pbmc.subset@meta.data[,c('condition', 'individual')])
    meta$cellCount <- sapply(meta$individual, function(id) sum(pbmc.subset$individual == id))
    meta <- meta[match(gsub('_.+', '', rownames(bulk)), meta$individual),]
    bulk <- cbind(class=meta$condition, individual=meta$individual, cellCount=meta$cellCount, bulk)
    cell <- gsub('/| |-', '.', cell)
    saveRDS(bulk, paste0('pseudobulk/', cell, '.HVG.RDS'))
}

# # Export the entire pseudobulked expression matrix, subsetted by X chromosome genes
# pbmc.subset <- subset(pbmc, features=rownames(chrX))
# bulk <- data.frame(t(AggregateExpression(pbmc.subset, group.by='individual', slot='counts')$RNA))
# meta <- unique(pbmc.subset@meta.data[,c('condition', 'individual')])
# meta <- meta[match(rownames(bulk), meta$individual),]
# bulk <- cbind(class=meta$condition, individual=meta$individual, bulk)
# saveRDS(bulk, 'psuedobulk/all.cells.chrX.RDS')

# # Export the entire pseudobulked expression matrix, subsetted by HVG
# pbmc.subset <- subset(pbmc, features = VariableFeatures(pbmc))
# bulk <- data.frame(t(AggregateExpression(pbmc.subset, group.by='individual', slot='counts')$RNA))
# meta <- unique(pbmc.subset@meta.data[,c('condition', 'individual')])
# meta <- meta[match(rownames(bulk), meta$individual),]
# bulk <- cbind(class=meta$condition, individual=meta$individual, bulk)
# saveRDS(bulk, 'psuedobulk/all.cells.HVG.RDS')

