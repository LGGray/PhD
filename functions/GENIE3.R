library(GENIE3)
library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')

# Load data, subset for celltype and condition then psuedobulk
pbmc <- readRDS('pbmc.female.RDS')

deg <- deg.list('differential.expression/edgeR', logfc=0.5)
names(deg) <- c(
  "CD16+ NK cells", "Classical monocytes", "DC1", "DC2", "MAIT cells", "Mast cells", "Memory B cells", 
  "Naive B cells", "NK cells", "Non-classical monocytes", "pDC", "Plasma cells", "Plasmablasts", "Regulatory T cells", 
  "Tcm/Naive cytotoxic T cells", "Tcm/Naive helper T cells", "Tem/Effector helper T cells", "Tem/Temra cytotoxic T cells", 
  "Tem/Trm cytotoxic T cells"
)

cell <- 'Memory B cells'
disease <- subset(pbmc, cellTypist == cell & condition == 'disease')
keep <- rowSums(disease@assays$decontXcounts@counts > 0) > ncol(disease) * 0.05
features <- names(keep[keep == T])
disease <- subset(disease, features=features)
# Psudobulking by summing counts
disease.expr <- AggregateExpression(disease, group.by='individual', slot='counts')$decontXcounts

TF <- read.delim('../../SCENIC/allTFs_hg38.txt', header=F)

GRTD <- clusterProfiler::read.gmt('../../gene.sets/c3.tft.gtrd.v7.5.1.symbols.gmt')

TF <- TF[TF$V1 %in% rownames(disease.expr),]
features <- rownames(chrX)[rownames(chrX) %in% rownames(disease.expr)]

# Run GENIE3
weightMat <- GENIE3(disease.expr, nCores=4, regulators=TF)
# Find regulatory links
linkList <- getLinkList(weightMat)
summary(linkList$weight)

rownames(weightMat)



weightMat.chrX <- GENIE3(disease.expr, nCores=4, regulators=TF, targets=features)
linkList.chrX <- getLinkList(weightMat.chrX)
summary(linkList.chrX$weight)

pdf('APR/GENIE3.chrX.heatmap.pdf')
Heatmap(weightMat.chrX, name='Weight', 
show_row_names = FALSE, show_column_names = FALSE)
dev.off()

