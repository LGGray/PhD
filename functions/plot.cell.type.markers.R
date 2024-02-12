library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

markers <- c('CD4', 'CD40LG', 'CD8A', 'CD8B', 'CCR7', 'IL7R',
'TNFRSF4', 'RTKN2', 'FOXP3', 'PRF1', 'GZMH', 'GZMB', 'GZMK', 
'KLRB1', 'GNLY', 'NKG7', 'FCGR3A', 'TCL1A', 'BANK1', 'MZB1', 
'FCRL5', 'CD19', 'MS4A1', 'CD79A', 'SOX4', 'CD14', 'FCGR3A', 
'C1QA', 'CLEC10A', 'CLEC9A', 'LILRA4')

load('/directflow/SCCGGroupShare/projects/lacgra/cellTypist/curated.markers.RData')
cellTypist <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/cellTypist/.txt', sep = '\t', header = T)



pbmc <- readRDS('pbmc.female.control-managed.RDS')

pseudobulk <- AverageExpression(pbmc, features = markers, slot='counts')$RNA
# pseudobulk <- scale(pseudobulk)

pseudobulk <- pseudobulk[rowSums(pseudobulk) > 0,]

pdf('marker_heatmap.pdf', width = 10, height = 10)
Heatmap(t(pseudobulk), name = 'Normalised expression', col = colorRamp2(c(0, 1), c('white', 'red')),
column_names_gp = gpar(fontsize = 12), row_names_gp = gpar(fontsize = 12))
dev.off()

# features <- rownames(pseudobulk)
# feature.celltype <- lapply(features, function(x) cellTypist[grep(x, cellTypist$Curated.markers), 'Low.hierarchy.cell.types'])
# names(feature.celltype) <- features

# pdf('marker.gene.dotplot.pdf', width = 20, height = 10)
# DotPlot(pbmc, features = curated.markers) + theme(axis.text.x = element_text(size=5))
# dev.off()