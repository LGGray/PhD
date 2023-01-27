library(Seurat)
library(SingleCellExperiment)
library(dplyr)
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')

args = commandArgs(trailingOnly=TRUE)

setwd(paste0('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/', args[1]))

if(dir.exists('cpdb') != T){dir.create('cpdb')}

deg <- edgeR.list('psuedobulk', logfc=0.5)
names(deg) <- gsub('.edgeR-LRT', '', names(deg))
names(deg) <- gsub('Tem_Trm', 'Tem/Trm', names(deg))
names(deg) <- gsub('Tem_Effector', 'Tem/Effector', names(deg))
deg.df <- bind_rows(deg, .id='celltype')
deg.df <- subset(deg.df, abs(logFC) > 0.5)
deg.df$celltype <- gsub('_', ' ', deg.df$celltype)

write.table(deg.df, 'cpdb/DEGs.tsv', sep='\t', row.names=F, quote = F)

pbmc <- readRDS('pbmc.female.RDS')
cells <- unique(deg.df$celltype)
pbmc <- subset(pbmc, cellTypist %in% cells)

writeMM(pbmc@assays$SCT@data, file = 'cpdb/matrix.mtx')
write(x = rownames(pbmc@assays$SCT@data), file = "cpdb/features.tsv")
write(x = colnames(pbmc@assays$SCT@data), file = "cpdb/barcodes.tsv")

pbmc@meta.data %>% 
  tibble::rownames_to_column(var='Cell') %>%
  select('Cell', 'cellTypist') %>%
  write.table(., 'cpdb/meta.tsv', sep='\t', quote=F, row.names = F)

