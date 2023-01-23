library(Seurat)
library(SeuratDisk)
library(sctransform)
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)

setwd('/directflow/SCCGGroupShare/projects/lacgra/lupus.Chun/')

args = commandArgs(trailingOnly=TRUE)
args <- as.numeric(args)

ancestry <- c('european', 'asian')

#Convert(paste0(ancestry[args], '.h5ad'), dest='h5seurat', overwrite=F)
pbmc <- LoadH5Seurat(paste0(ancestry[args], '.h5seurat'), 
                       meta.data=F, assays='RNA', misc=F)

# Add metadata
metadata <- read.csv(paste0(ancestry[args],'.metadata.csv'))
rownames(metadata) <- metadata$index
metadata <- metadata[,-1]
pbmc@meta.data <- metadata

expr <- GetAssayData(pbmc, assay='RNA')

rm(pbmc)

converted <- select(EnsDb.Hsapiens.v86, # database
                    keys = rownames(expr),  # data to use for retrieval
                    column = "SYMBOL", # information to retreive for given data
                    keytype = "GENEID")

rownames(expr) <- converted$SYMBOL

pbmc <- CreateSeuratObject(expr)

pbmc@meta.data <- cbind(pbmc@meta.data, metadata)

# Transform data
pbmc <- SCTransform(pbmc)

# Read in PBMC reference dataset
reference <- LoadH5Seurat("~/azimuth.reference/pbmc_multimodal.h5seurat")
# Find anchors between reference and query
anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)
# Transfer cell type labels and protein data from reference to query
pbmc <- MapQuery(
  anchorset = anchors,
  query = pbmc,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

Idents(pbmc) <- 'predicted.celltype.l2'

saveRDS(pbmc, paste0(ancestry[args], '.RDS'))

# subset female
pbmc.female <- subset(pbmc, sex == 'female')
saveRDS(pbmc, paste0(ancestry[args], '.female.RDS'))

# load('~/datasets/integrated/feature.selection/randomForest/CD16_Mono.RData')
# xcape <- read.delim('~/datasets/integrated/feature.selection/randomForest/CD16_Mono.txt')
# features <- xcape$gene
# 
# pbmc.subset <- subset(pbmc, cell_type %in% 'non-classical monocyte')
# expr <- t(as.matrix(GetAssayData(pbmc.subset, assay='RNA')))
# Convert ENSG to Symbol
# converted <- select(EnsDb.Hsapiens.v86, # database
#                     keys = rownames(expr),  # data to use for retrieval
#                     column = "SYMBOL", # information to retreive for given data
#                     keytype = "GENEID")
# 
# rownames(expr) <- converted$SYMBOL
# expr <- data.frame(expr[,which(colnames(expr) %in% features)])
# expr$condition <- factor(pbmc.subset$condition, levels = c('Control', 'Disease'))
