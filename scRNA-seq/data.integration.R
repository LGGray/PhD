library(Seurat)
library(ggplot2)
source('R_code/colour.dictionary.R')


setwd('datasets')

# Read in objects
RA <- readRDS('SDY998/pbmc.female.RDS')
SLE <- readRDS('SDY997/pbmc.lymphocytes.female.RDS')
UC <- readRDS('GSE125527/pbmc.female.RDS')
pSS <- readRDS('GSE157278/pbmc.RDS')
MS <- readRDS('GSE19')

UC$disease_assignment <- gsub('diseased', 'UC', UC$disease_assignment)
UC$disease_assignment <- gsub('healthy', 'Control', UC$disease_assignment)
colnames(UC@meta.data)[7] <- 'disease'

colnames(pSS@meta.data)[5] <- 'disease'
pSS$disease <- gsub('HC', 'Control', pSS$disease)

colnames(MS@meta.data)[15] <- 'disease'
MS$disease <- gsub('HC', 'Control', MS$disease)

object.list <- list(RA, SLE, UC, pSS, MS)
names(object.list) <- c('RA', 'SLE', 'UC', 'pSS', 'MS')

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = 3000)
object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = features)

# Identify anchors 
immune.anchors <- FindIntegrationAnchors(object.list = object.list, 
                                         normalization.method = "SCT",
                                         anchor.features = features)

# Integrated data assay
immune.combined <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

immune.combined <- RunPCA(immune.combined, verbose = FALSE)
ElbowPlot(immune.combined)
dev.off()
immune.combined <- FindNeighbors(immune.combined, dims=1:15)
immune.combined <- FindClusters(immune.combined)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims=1:15)

Idents(immune.combined) <- 'predicted.celltype.l2'

pdf('integrated/UMAP.integrated.pdf')
cells <- gsub(' ', '_', levels(immune.combined))
DimPlot(immune.combined, cols = as.vector(unlist(cell.colourdict[cells])), label=T)
dev.off()

saveRDS(immune.combined, 'integrated/immune.combined.RDS')
