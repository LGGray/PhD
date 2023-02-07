library(Seurat)
library(Matrix)
library(dplyr)
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')

args = commandArgs(trailingOnly=TRUE)

setwd(paste0('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/', args[1]))

if(dir.exists('cpdb') != T){dir.create('cpdb')}

deg <- edgeR.list('psuedobulk', logfc=0.5)
deg <- lapply(deg, function(x) subset(x, logFC > 0.5))
names(deg) <- gsub('.edgeR-LRT', '', names(deg))
deg.df <- bind_rows(deg, .id='cellTypist')

write.table(deg.df, 'cpdb/DEGs.tsv', sep='\t', row.names=F, quote = F)

pbmc <- readRDS('pbmc.female.RDS')
pbmc@meta.data$cellTypist <- gsub('/|-| ', '_', pbmc@meta.data$cellTypist)
Idents(pbmc) <- pbmc@meta.data$cellTypist
cells <- unique(deg.df$cellTypist)
pbmc <- subset(pbmc, cellTypist %in% cells)

writeMM(pbmc@assays$SCT@data, file = 'cpdb/matrix.mtx')
write(x = rownames(pbmc@assays$SCT@data), file = "cpdb/features.tsv")
write(x = colnames(pbmc@assays$SCT@data), file = "cpdb/barcodes.tsv")

pbmc@meta.data$Cell = rownames(pbmc@meta.data)
df = pbmc@meta.data[, c('Cell', 'cellTypist')]
write.table(df, file ='cpdb/meta.tsv', sep = '\t', quote = F, row.names = F)