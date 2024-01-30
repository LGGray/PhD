library(Seurat)
library(ggplot2)
library(reshape2)

pbmc <- readRDS('pbmc.female.control-managed.RDS')

markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, 
return.thresh = 0.8, test.use = 'roc', slot='scale.data')

# save(markers, file='celltype.markers.RData')
load('celltype.markers.RData')

marker.genes <- unique(subset(markers, myAUC > 0.9)$gene)
pseudobulked <- AverageExpression(pbmc, features=marker.genes, group.by='cellTypist', slot='counts')$RNA
pseudobulked <- pseudobulked[rowSums(pseudobulked) > 0,]
marker.genes <- rownames(pseudobulked)
pseudobulked.long <- melt(pseudobulked, id.vars='cellTypist', value.name='expression')

pdf('celltype.markers.pdf', width=10, height=10)
ggplot(pseudobulked.long, aes(x=Var1, y=Var2, colour=expression)) + 
    geom_point() +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(x='Cell type', y='Average expression', fill='Marker gene')
dev.off()

pdf('celltype.markers.dotplot.pdf', width=20, height=10)
DotPlot(pbmc, features = marker.genes) + RotatedAxis()
dev.off()

pdf('celltype.markers.heatmap.pdf', width=10, height=10)
DoHeatmap(pbmc, features = marker.genes, group.by='cellTypist', slot='scale.data')
dev.off()