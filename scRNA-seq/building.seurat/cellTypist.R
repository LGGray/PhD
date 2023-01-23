library(Seurat)

args <- commandArgs()

setwd(paste0('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets', args))

pbmc <- readRDS('pbmc.RDS')

labels <- read.csv('cellTypist/predicted_labels.csv')
pbmc@meta.data <- cbind(pbmc@meta.data, cellTypist=labels$majority_voting)

unique(paste(pbmc$cellTypist, pbmc$seurat_clusters))

Idents(pbmc) <- 'cellTypist'

pdf('DimPlot.cellTypist.pdf')
DimPlot(pbmc, label = TRUE, reduction='ref.umap', repel=T) + NoLegend()
dev.off()

pbmc.markers <- FindAllMarkers(pbmc, only.pos=T, min.pct=0.25, logfc.threshold = 0.25)

write.table(pbmc.markers, 'cellTypist.markers.txt')

