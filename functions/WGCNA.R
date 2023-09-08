library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)
library(WGCNA)
library(dynamicTreeCut)

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')

# Load data and psuedobulk
pbmc <- readRDS('pbmc.female.control-managed.RDS')

cell = 'Tem/Temra cytotoxic T cells'
pbmc.cell <- subset(pbmc, cellTypist == cell)

# Keep genes with expression in 5% of cells
keep <- rowSums(pbmc.cell@assays$RNA@counts > 0) > ncol(pbmc.cell) * 0.05
features <- names(keep[keep == T])
pbmc.cell <- subset(pbmc.cell, features=features)

# Psudobulking by summing counts
expr <- AggregateExpression(pbmc.cell, group.by='individual', slot='counts')$RNA
expr <- data.frame(t(expr))

gsg = goodSamplesGenes(expr, verbose = 3)
if (!gsg$allOK)
{
# Optionally, print the gene and sample names that were removed:
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(expr)[!gsg$goodGenes], collapse = ", ")));
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(expr)[!gsg$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
expr = expr[gsg$goodSamples, gsg$goodGenes]
}

# WGCNA analysis
options(stringsAsFactors = FALSE)
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=200, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(expr, powerVector = powers, verbose = 5, networkType = "signed")
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
pdf('scale_free_topology.pdf')
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
dev.off()
# Mean connectivity as a function of the soft-thresholding power
pdf('mean_connectivity.pdf')
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# Calculate adjecency matrix
softPower = 200
adjacency = adjacency(expr, power = softPower)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
pdf('gene_dendrogram.pdf')
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04)
dev.off()

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
pdf('module assignments.pdf')
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(expr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
sizeGrWindow(7, 6)
pdf('module_eigengenes.pdf')
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;