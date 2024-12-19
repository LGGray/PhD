library(Seurat)
library(edgeR)

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/Seurat2PB.R')

# Load X chromosome genes
chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name

# Load SLE DisGenNet genes
SLE <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/SLE.tsv')
SLE <- unique(SLE$Gene)

# Function to calculate gene variability
calculate_variability <- function(expression_data) {
  variability <- apply(expression_data, 2, var)
  return(variability)
}

# Read in Seurat object
pbmc <- readRDS('pbmc_female.control_managed.RDS')

###### Pseudobulk analysis ######

# Export psuedobulked expression matrix, subsetted by X chromosome and autosomal genes for each cell type
for(cell in levels(pbmc)){
    pbmc.subset <- subset(pbmc, cell_type_detailed == cell)
    # Pseudobulk by summed expression
    bulk <- AggregateExpression(pbmc.subset, slot='counts', group.by='ind_cov')$COMBAT_LogNorm
    # Remove lowly expressed genes
    non_zero_genes <- rowSums(bulk) > 0
    bulk <- bulk[non_zero_genes, ]
    # Divide each gene by total number of counts across all genes
    bulk <- apply(bulk, 2, function(x) x / sum(x))

    # Subset to X chromosome genes
    bulk.chrX <- data.frame(t(bulk[rownames(bulk) %in% chrX,]))
    # Subset to autosomal genes
    bulk.autosome <- data.frame(t(bulk[!rownames(bulk) %in% chrX,]))
    # Calculate gene variability
    chrX.variability <- calculate_variability(bulk.chrX)
    autosomal.variability <- calculate_variability(bulk.autosome)
    # Match autosomal genes to X chromosome genes by variability
    matched_autosomal_genes <- vector()
    for (var in chrX.variability) {
    closest_gene <- names(which.min(abs(autosomal.variability - var)))
    matched_autosomal_genes <- c(matched_autosomal_genes, closest_gene)
    autosomal.variability <- autosomal.variability[!names(autosomal.variability) %in% closest_gene]
    }
    matched_autosomal_genes <- unique(matched_autosomal_genes)
    bulk.autosome <- bulk.autosome[,matched_autosomal_genes]

    # Create metadata to add to pseudobulk matrix
    meta <- unique(pbmc.subset@meta.data[,c('disease', 'ind_cov', 'age', 'ancestry')])
    meta$cellCount <- sapply(meta$ind_cov, function(id) sum(pbmc.subset$ind_cov == id))
    meta <- meta[match(rownames(bulk.chrX), meta$ind_cov),]

    # Add metadata to each pseudobulk matrix
    bulk.chrX <- cbind(class=meta$disease, individual=meta$ind_cov, cellCount=meta$cellCount, age=meta$age,
    ancestry=meta$ancestry, bulk.chrX)
    bulk.autosome <- cbind(class=meta$disease, individual=meta$ind_cov, cellCount=meta$cellCount, age=meta$age,
    ancestry=meta$ancestry, bulk.autosome)
    cell <- gsub("\\+|-| ", "_", cell)
    saveRDS(bulk.chrX, paste0('pseudobulk/', cell, '.chrX.RDS'))
    saveRDS(bulk.autosome, paste0('pseudobulk/', cell, '.autosome.RDS'))
}

# Export the pseudobulked expression matrix, subsetted by HVG for each cell type
for(cell in levels(pbmc)){
    pbmc.subset <- subset(pbmc, cell_type_detailed == cell)
    # Pseudobulk by summed expression
    bulk <- AggregateExpression(pbmc.subset, slot='counts', group.by='ind_cov')$COMBAT_LogNorm
    # Remove lowly expressed genes
    non_zero_genes <- rowSums(bulk) > 0
    bulk <- bulk[non_zero_genes, ]
    # Divide each gene by total number of counts across all genes
    bulk <- apply(bulk, 2, function(x) x / sum(x))

    # Calulate variation for each gene and select top N genes
    # N is the nuber of X chromosome genes
    variation <- apply(bulk, 1, var)
    variation <- sort(variation, decreasing=TRUE)
    hvg.all <- names(variation[1:sum(chrX %in% rownames(bulk))])
    # Remove X chromosome genes and select top N autosomal genes
    variation <- variation[!names(variation) %in% chrX]
    hvg.autosome <- names(variation[1:sum(chrX %in% rownames(bulk))])

    # Subset data
    bulk.HVG <- data.frame(t(bulk[rownames(bulk) %in% hvg.all,]))
    bulk.HVG.autosome <- data.frame(t(bulk[rownames(bulk) %in% hvg.autosome,]))

    # Create metadata to add to pseudobulk matrix
    meta <- unique(pbmc.subset@meta.data[,c('disease', 'ind_cov', 'age', 'ancestry')])
    meta$cellCount <- sapply(meta$ind_cov, function(id) sum(pbmc.subset$ind_cov == id))
    meta <- meta[match(rownames(bulk.HVG), meta$ind_cov),]

    # Add metadata to each pseudobulk matrix
    bulk.HVG <- cbind(class=meta$disease, individual=meta$ind_cov, cellCount=meta$cellCount, age=meta$age,
    ancestry=meta$ancestry, bulk.HVG)
    bulk.HVG.autosome <- cbind(class=meta$disease, individual=meta$ind_cov, cellCount=meta$cellCount, age=meta$age,
    ancestry=meta$ancestry, bulk.HVG.autosome)
    cell <- gsub("\\+|-| ", "_", cell)
    saveRDS(bulk.HVG, paste0('pseudobulk/', cell, '.HVG.RDS'))
    saveRDS(bulk.HVG.autosome, paste0('pseudobulk/', cell, '.HVG.autosome.RDS')) 
}

# Export the pseudobulked expression matrix, subsetted by SLE DisGenNet
for(cell in levels(pbmc)){
    pbmc.subset <- subset(pbmc, cell_type_detailed == cell)
    # Pseudobulk by summed expression
    bulk <- AggregateExpression(pbmc.subset, slot='counts', group.by='ind_cov')$COMBAT_LogNorm
    # Remove lowly expressed genes
    non_zero_genes <- rowSums(bulk) > 0
    bulk <- bulk[non_zero_genes, ]
    # Divide each gene by total number of counts across all genes
    bulk <- apply(bulk, 2, function(x) x / sum(x))

    # Subset for SLE genes
    bulk.SLE <- data.frame(t(bulk[rownames(bulk) %in% SLE,]))

    # Create metadata to add to pseudobulk matrix
    meta <- unique(pbmc.subset@meta.data[,c('disease', 'ind_cov', 'age', 'ancestry')])
    meta$cellCount <- sapply(meta$ind_cov, function(id) sum(pbmc.subset$ind_cov == id))
    meta <- meta[match(rownames(bulk.SLE), meta$ind_cov),]

    # Add metadata to each pseudobulk matrix
    bulk.SLE <- cbind(class=meta$disease, individual=meta$ind_cov, cellCount=meta$cellCount, age=meta$age,
    ancestry=meta$ancestry, bulk.SLE)
    cell <- gsub("\\+|-| ", "_", cell)
    saveRDS(bulk.SLE, paste0('pseudobulk/', cell, '.SLE.RDS'))
}

