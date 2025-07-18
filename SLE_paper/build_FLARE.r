library(reticulate)
library(Matrix)
library(Seurat)
library(scico)
library(dplyr)

scipy_sparse <- import("scipy.sparse")
mtx <- scipy_sparse$load_npz("FLARE/exp_counts_female_flare.npz")
mtx <- t(mtx)

# Read in metadata
metadata <- read.delim('FLARE/metadata.tsv', header=T, row.names=1)
metadata$development_stage <- gsub('-year-old stage', '', metadata$development_stage)
metadata$development_stage <- as.numeric(metadata$development_stage)
metadata$disease <- ifelse(metadata$disease == 'systemic lupus erythematosus', 'disease', 'control')
colnames(metadata)[20] <- 'batch'
colnames(metadata)[30] <- 'ancestry'
colnames(metadata)[31] <- 'age'


features <- read.delim('FLARE/features.tsv', header=F)
colnames(features) <- c('gene_id', 'feature_name', 'feature_biotype')
# barcodes <- read.delim('barcodes.tsv.gz', header=F)
rownames(mtx) <- features$feature_name
colnames(mtx) <- rownames(metadata)

# Create a Seurat object without normalizing again
pbmc <- CreateSeuratObject(counts = mtx, assay = "COMBAT_LogNorm", min.cells = 0, min.features = 0)

### Add the raw counts matrix to the Seurat object
mtx <- scipy_sparse$load_npz("FLARE/raw_counts_female_flare.npz")
mtx <- t(mtx)
rownames(mtx) <- features$feature_name
colnames(mtx) <- rownames(metadata)

pbmc[["RNA"]] <- CreateAssayObject(counts = mtx, min.cells = 0, min.features = 0)

# Set the data slot to log-normalized values
DefaultAssay(pbmc) <- "RNA"

# Add metadata to Seurat object
pbmc@meta.data <- cbind(pbmc@meta.data, metadata)

# Cell type clustering
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunPCA(pbmc)

# Find the number of PCs to use
# Plot the elbow plot
pdf('FLARE/seurat.pca.pdf')
ElbowPlot(pbmc, ndims = 50)
dev.off()

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:11)

# Create detailed cell type column
pbmc@meta.data$cell_type_detailed <- case_when(
  # CD4 T cell subsets
  pbmc$cell_type == "CD4-positive, alpha-beta T cell" & pbmc$ct_cov == "T4_naive" ~ "CD4 T cell_Naive",
  pbmc$cell_type == "CD4-positive, alpha-beta T cell" & pbmc$ct_cov == "T4_em"    ~ "CD4 T cell_Effector Memory",
  pbmc$cell_type == "CD4-positive, alpha-beta T cell" & pbmc$ct_cov == "T4_reg"   ~ "CD4 T cell_Treg",

  # CD8 T cell subsets
  pbmc$cell_type == "CD8-positive, alpha-beta T cell" & pbmc$ct_cov == "T8_naive" ~ "CD8 T cell_Naive",
  pbmc$cell_type == "CD8-positive, alpha-beta T cell" & pbmc$ct_cov == "CytoT_GZMH+" ~ "CD8 T cell_Cytotoxic GZMH+",
  pbmc$cell_type == "CD8-positive, alpha-beta T cell" & pbmc$ct_cov == "CytoT_GZMK+" ~ "CD8 T cell_Cytotoxic GZMK+",
  pbmc$cell_type == "CD8-positive, alpha-beta T cell" & pbmc$ct_cov == "T_mait"   ~ "CD8 T cell_MAIT",

  # NK cell subsets
  pbmc$cell_type == "natural killer cell" & pbmc$ct_cov == "NK_bright" ~ "NK cell_Bright",
  pbmc$cell_type == "natural killer cell" & pbmc$ct_cov == "NK_dim"    ~ "NK cell_Dim",

  # B cell subsets
  pbmc$cell_type == "B cell" & pbmc$ct_cov == "B_naive"    ~ "B cell_Naive",
  pbmc$cell_type == "B cell" & pbmc$ct_cov == "B_mem"      ~ "B cell_Memory",
  pbmc$cell_type == "B cell" & pbmc$ct_cov == "B_plasma"   ~ "B cell_Plasma",
  pbmc$cell_type == "B cell" & pbmc$ct_cov == "B_atypical" ~ "B cell_Atypical",

  # Progenitor cells
  pbmc$cell_type == "progenitor cell" ~ "Progenitor cell",

  # Dendritic cells
  pbmc$cell_type == "conventional dendritic cell" ~ "Conventional DC",
  pbmc$cell_type == "plasmacytoid dendritic cell" ~ "Plasmacytoid DC",

  # Monocytes
  pbmc$cell_type == "classical monocyte"       ~ "Classical monocyte",
  pbmc$cell_type == "non-classical monocyte"   ~ "Non-classical monocyte",

  # Fallback for lymphocyte and plasmablast without ct_cov
  pbmc$cell_type == "lymphocyte"               ~ "Lymphocyte",
  pbmc$cell_type == "plasmablast"              ~ "Plasmablast",

  # Default to cell_type for anything else
  TRUE ~ pbmc$cell_type
)

# Cell type colours
cell_type_colours <- c(
  # T cells - shades of blue and orange
  "CD4 T cell_Naive" = "#377EB8",
  "CD4 T cell_Effector Memory" = "#6AAED6",
  "CD4 T cell_Treg" = "#1F78B4",
  "CD8 T cell_Naive" = "#FF7F00",
  "CD8 T cell_Cytotoxic GZMH+" = "#FDBF6F",
  "CD8 T cell_Cytotoxic GZMK+" = "#FE9929",
  "CD8 T cell_MAIT" = "#E6550D",

  # NK cells - shades of red
  "NK cell_Bright" = "#E31A1C",
  "NK cell_Dim" = "#FC9272",
  "natural killer cell" = "#FB6A4A",

  # B cells - shades of purple
  "B cell_Naive" = "#CAB2D6",
  "B cell_Memory" = "#9E9AC8",
  "B cell_Plasma" = "#6A3D9A",
  "B cell_Atypical" = "#BC80BD",

  # Monocytes - shades of brown
  "Classical monocyte" = "#8C510A",
  "Non-classical monocyte" = "#D8B365",

  # Dendritic cells - shades of green
  "Conventional DC" = "#33A02C",
  "Plasmacytoid DC" = "#B2DF8A",

  # Other - distinct colors
  "Lymphocyte" = "#A6CEE3",
  "Plasmablast" = "#FB9A99",
  "Progenitor cell" = "#E6AB02"
)

# Plot UMAP
Idents(pbmc) <- 'cell_type_detailed'

pdf('FLARE/seurat.umap.pdf', width = 10, height = 10)
DimPlot(pbmc, cols = cell_type_colours, label = TRUE, 
label.size = 5, repel = TRUE)
dev.off()

# Save Seurat object
saveRDS(pbmc, 'FLARE/pbmc_female.control_FLARE.RDS')

###### Export Pseudobulk #####
library(edgeR)
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/Seurat2PB.R')

# Load X chromosome genes
chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name

# Load SLE DisGenNet genes
SLE <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/SLE.tsv')
SLE <- unique(SLE$Gene)

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

    # Create metadata to add to pseudobulk matrix
    meta <- unique(pbmc.subset@meta.data[,c('disease', 'ind_cov', 'age', 'ancestry')])
    meta$cellCount <- sapply(meta$ind_cov, function(id) sum(pbmc.subset$ind_cov == id))
    meta <- meta[match(rownames(bulk.chrX), meta$ind_cov),]

    # Add metadata to each pseudobulk matrix
    bulk.chrX <- cbind(class=meta$disease, individual=meta$ind_cov, cellCount=meta$cellCount, age=meta$age,
    ancestry=meta$ancestry, bulk.chrX)

    # Save pseudobulk matrix
    cell <- gsub("\\+|-| ", "_", cell)
    saveRDS(bulk.chrX, paste0('FLARE/pseudobulk/', cell, '.chrX.RDS'))
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
    saveRDS(bulk.SLE, paste0('FLARE/pseudobulk/', cell, '.SLE.RDS'))
}