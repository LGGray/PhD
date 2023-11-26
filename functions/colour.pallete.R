# Script to generate a colour pallete for all cell types studied
# We create a unique list of cell types from the ouput of propeller and assign a colour to the cell type

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')

pSS <- read.delim('pSS_GSE157278/propeller.asin.txt', sep=' ')
SLE <- read.delim('lupus_Chun/propellor.asin.age.condition.abundance.txt')
UC <- read.delim('UC_GSE125527/propeller.asin.txt', sep=' ')
CO <- readRDS('CD_Kong/colon/propellor.asin.RDS')
TI <- readRDS('CD_Kong/TI/propellor.asin.RDS')

celltypes <- unique(c(rownames(pSS), rownames(SLE), rownames(UC), rownames(CO), rownames(TI)))
# Reorder cell types based on lineage
celltypes <- c("B cells", "Naive B cells", "Memory B cells", "Age-associated B cells", "Follicular B cells", 
"Germinal center B cells", "Proliferative germinal center B cells", "Plasmablasts", "Plasma cells", "Double-positive thymocytes", 
"Tcm/Naive cytotoxic T cells", "Tem/Trm cytotoxic T cells", "Tem/Temra cytotoxic T cells", "Trm cytotoxic T cells", 
"Tcm/Naive helper T cells", "Tem/Effector helper T cells", "Type 1 helper T cells", "Type 17 helper T cells", "Regulatory T cells", 
"Follicular helper T cells", "Cycling T cells", "MAIT cells", "gamma-delta T cells", "CRTAM+ gamma-delta T cells", "NK cells", 
"CD16+ NK cells", "CD16- NK cells", "ILC", "ILC3", "Monocytes", "Classical monocytes", "Non-classical monocytes", "Macrophages", 
"Intestinal macrophages", "Intermediate macrophages", "Erythrophagocytic macrophages", "DC1", "DC2", "pDC", "Migratory DCs", 
"Mast cells", "Myelocytes", "HSC/MPP")

# Create a colourblind safe pallete for each cell type
library(viridis)

colours <- viridis(length(celltypes))
names(colours) <- celltypes

# Save the colour pallete
save(colours, file='celltype.colours.RData')

load('../../celltype.colours.RData')
# DimPlot
pdf('DimPlot.cellTypist.female.pdf')
DimPlot(pbmc, reduction='umap', group.by='cellTypist',label.color = "black", label=TRUE, repel=TRUE, pt.size=0.5, cols=colours) + NoLegend()
dev.off()