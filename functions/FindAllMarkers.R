library(Seurat)

pbmc <- readRDS('pbmc.female.control-managed.RDS')

markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, 
return.thresh = 0.8, test.use = 'roc', slot='scale.data')

save(markers, file='celltype.markers.RDS')