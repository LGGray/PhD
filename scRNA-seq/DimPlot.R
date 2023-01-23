library(Seurat)
library(RColorBrewer)
library(ggplot2)
source('R_code/colour.dictionary.R')
load('XCI/escapees.Rdata')
load('XCI/chrX.Rdata')


setwd('datasets')

# Read in objects
RA <- readRDS('SDY998/pbmc.female.RDS')
SLE <- readRDS('SDY997/pbmc.lymphocytes.female.RDS')
UC <- readRDS('GSE125527/pbmc.female.RDS')
pSS <- readRDS('GSE157278/pbmc.RDS')
MS <- readRDS('GSE193770/pbmc.female.RDS')

# Colour by cell type
cells <- gsub(' ', '_', levels(MS))
pdf('GSE193770/UMAP.celltype.pdf')
DimPlot(MS, cols = as.vector(unlist(cell.colourdict[cells])), label=T)
dev.off()

# Colour by condition
pdf('GSE193770/UMAP.condition.pdf')
DimPlot(MS, group.by = 'condition')
dev.off()
