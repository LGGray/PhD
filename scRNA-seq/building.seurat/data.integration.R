library(Seurat)
library(ggplot2)

setwd('datasets')

# Read in objects
RA <- readRDS('SDY998/pbmc.female.RDS')
SLE <- readRDS('SDY997/pbmc.lymphocytes.female.RDS')
UC <- readRDS('GSE125527/pbmc.female.RDS')
pSS <- readRDS('GSE157278/pbmc.RDS')

UC$disease_assignment <- gsub('diseased', 'UC', UC$disease_assignment)

experimental.condition <- data.frame
experimental.condition <- cbind(experimental.condition, RA$disease, 
                                SLE$disease, UC$disease_assignment, pSS$condition)

object.list <- list(RA, SLE, UC, pSS)'disease'
names(object.list) <- c('RA', 'SLE', 'UC', 'pSS')

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
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)

# Pull out condition information
condition <- metadata[,c('disease', 'disease_assignment', 'condition')]

a <- condition[!is.na(condition[,1]),1]
b <- condition[!is.na(condition[,2]),2]
c <- condition[!is.na(condition[,3]),3]

condition$integrated.condition <- c(a,b,c)

immune.combined@meta.data$integrated.condition <- condition$integrated.condition

DimPlot(immune.combined, reduction = "umap", group.by='predicted.celltype.l2', label=T)
dev.off()

saveRDS(immune.combined, 'integrated/immune.combined.RDS')
