library(reticulate)
library(Matrix)
library(Seurat)
library(scico)

scipy_sparse <- import("scipy.sparse")
mtx <- scipy_sparse$load_npz("exp_counts_female_managed.npz")
mtx <- t(mtx)
features <- read.delim('features.tsv.gz', header=F)
colnames(features) <- c('gene_id', 'feature_name', 'feature_biotype')
barcodes <- read.delim('barcodes.tsv.gz', header=F)
rownames(mtx) <- features$feature_name
colnames(mtx) <- barcodes$V1

# Create a Seurat object without normalizing again
pbmc <- CreateSeuratObject(counts = mtx, assay = "COMBAT_LogNorm", min.cells = 0, min.features = 0)

# Set the data slot to log-normalized values
DefaultAssay(pbmc) <- "COMBAT_LogNorm"

# Read in metadata
metadata <- read.delim('metadata.tsv', header=T, row.names=1)
metadata$development_stage <- gsub('-year-old stage', '', metadata$development_stage)
metadata$development_stage <- as.numeric(metadata$development_stage)
metadata$disease <- ifelse(metadata$disease == 'systemic lupus erythematosus', 'disease', 'control')
colnames(metadata)[20] <- 'batch'
colnames(metadata)[30] <- 'ancestry'
colnames(metadata)[31] <- 'age'

# Add metadata to Seurat object
pbmc@meta.data <- cbind(pbmc@meta.data, metadata)

# Cell type clustering
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunPCA(pbmc)

# Find the number of PCs to use
# Plot the elbow plot
pdf('seurat.pca.pdf')
ElbowPlot(pbmc, ndims = 50)
dev.off()

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:11)

# Plot UMAP
Idents(pbmc) <- 'cell_type'
cell_type_colours <- scico(length(levels(pbmc)), palette = 'roma')
names(cell_type_colours) <- levels(pbmc$cell_type)
pdf('figures/seurat.umap.pdf')
DimPlot(pbmc, cols = cell_type_colours, label = FALSE)
dev.off()

# Save Seurat object
saveRDS(pbmc, 'pbmc_female.control_managed.RDS')
