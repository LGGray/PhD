# Script for building Seurat object from UC_GSE182270 data

library(Seurat)
library(magrittr)

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/UC_JamesLab')

pbmc <- readRDS('pbmc.unlabelled.RDS')

labels <- read.csv('cellTypist/predicted_labels.csv')
pbmc@meta.data <- cbind(pbmc@meta.data, cellTypist=labels$majority_voting)

Idents(pbmc) <- 'cellTypist'

pdf('DimPlot.cellTypist.all.pdf')
DimPlot(pbmc, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()

# # Determine sex of individuals by psuedobulked expression of chrY and XIST
# # pseduobulk expression matrix
# exp <- AverageExpression(pbmc, assays='SCT', features=c('XIST', rownames(chrY)), group.by='individual')$SCT
# exp <- scale(exp)

# # First we infer sex based on expression of female specific XIST gene
# XIST.expression <- exp[grep('XIST', rownames(exp)),]
# print(XIST.expression)

# # Perform hierarchical clustering to identify groups
# dissimilarity <- dist(data.frame(t(exp)), method='euclidean')
# cluster <- hclust(dissimilarity, method = 'centroid')

# # Plot dendrogram
# pdf('sex.dendrogram.pdf')
# plot(cluster)
# dev.off()

# # K-means clustering on the hclust data
# cluster.result <- cutree(cluster, k=2)
# # Assign sex based on dendrogram
# sex.list <- ifelse(cluster.result == 2, 'M', 'F')
# # Add sex to metadata
# pbmc$sex <- sex.list[pbmc$individual]

# Output all cells
DefaultAssay(pbmc) <- 'SCT'
saveRDS(pbmc, 'pbmc.RDS')

# Subset for females and output 
pbmc.female <- subset(pbmc, sex == 'female')

pdf('DimPlot.female.pdf')
DimPlot(pbmc.female, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()

saveRDS(pbmc.female, 'pbmc.female.RDS')

# Subset for males and output 
pbmc.male <- subset(pbmc, sex == 'male')

pdf('DimPlot.male.pdf')
DimPlot(pbmc.male, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()

saveRDS(pbmc.male, 'pbmc.male.RDS')