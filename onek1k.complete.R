library(Seurat)

pbmc <- readRDS("/directflow/SCCGGroupShare/projects/sarba2/data/onek1k/samples/cell_type.RDS")
load("/directflow/SCCGGroupShare/projects/sarba2/data/onek1k/samples/metadata.Rdata")
pbmc@meta.data = metadata

Idents(pbmc) <- 'predicted.celltype.l3'
pbmc$sex <- ifelse(pbmc$sex == 1, 'M', 'F')
print(levels(pbmc))
print(length(levels(pbmc)))
saveRDS(pbmc, '/directflow/SCCGGroupShare/projects/lacgra/seurat.object/onek1k.RDS')
