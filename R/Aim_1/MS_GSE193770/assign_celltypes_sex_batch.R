# Script for building Seurat object from SDY998 data

library(Seurat)
library(magrittr)

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrY.Rdata')

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/MS_GSE193770')

pbmc <- readRDS('pbmc.unlabelled.RDS')

labels <- read.csv('cellTypist/predicted_labels.csv')
pbmc@meta.data <- cbind(pbmc@meta.data, cellTypist=labels$majority_voting)

Idents(pbmc) <- 'cellTypist'

# Load in SVA output
load("svaseq.RData")
pbmc$batch_1 <- svseq$sv[,1]
pbmc$batch_2 <- svseq$sv[,2]

pdf('DimPlot.cellTypist.all.pdf')
DimPlot(pbmc, label = TRUE, reduction='umap', repel=T)
dev.off()

# Determine sex of individuals by psuedobulked expression of chrY and XIST
# pseduobulk expression matrix
exp <- AverageExpression(pbmc, assays='RNA', slot='counts', features=c('XIST', rownames(chrY)), group.by='individual')$RNA
exp <- scale(exp)

# First we infer sex based on expression of female specific XIST gene
XIST.expression <- exp[grep('XIST', rownames(exp)),]

# Perform hierarchical clustering to identify groups
dissimilarity <- dist(data.frame(t(exp)), method='euclidean')
cluster <- hclust(dissimilarity, method = 'centroid')

# Plot dendrogram
pdf('sex.dendrogram.pdf')
plot(cluster)
dev.off()

# K-means clustering on the hclust data
cluster.result <- cutree(cluster, k=2)
# Check the differencee between in expression of XIST between the two clusters
xist.1 <- mean(XIST.expression[names(which(cluster.result==1))])
xist.2 <- mean(XIST.expression[names(which(cluster.result==2))])
# Assign sex based on dendrogram
if(xist.1 > xist.2){
    sex.list <- ifelse(cluster.result == 1, 'F', 'M')
} else{
    sex.list <- ifelse(cluster.result == 1, 'M', 'F')
}
# Add sex to metadata
pbmc$sex <- sex.list[pbmc$individual]

# Output all cells
saveRDS(pbmc, 'pbmc.RDS')

# Subset for females and output 
pbmc.female <- subset(pbmc, sex == 'F')

pdf('DimPlot.female.pdf')
DimPlot(pbmc.female, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()
saveRDS(pbmc.female, 'pbmc.female.RDS')

# Subset for males and output 
pbmc.male <- subset(pbmc, sex == 'M')

pdf('DimPlot.male.pdf')
DimPlot(pbmc.male, label = TRUE, reduction='umap', repel=T) + NoLegend()
dev.off()
saveRDS(pbmc.male, 'pbmc.male.RDS')