library("UpSetR")

setwd('~/external/ClusterHome/datasets')
source('R_code/functions/edgeR.list.R')
load('XCI/chrX.Rdata')
load('XCI/escapees.Rdata')

degs <- edgeR.list('datasets/GSE157278/edgeR-LRT/', logfc=0.5)

# Filter genes for X chromosome
Xchr <- lapply(degs, function(x){
  subset(x, gene %in% rownames(chrX))$gene
})
upset(fromList(Xchr), nsets = length(names(Xchr)), nintersects = NA)
# Filter genes for XCI escape
xcape <- lapply(degs, function(x){
  subset(x, gene %in% rownames(escape))$gene
})
upset(fromList(xcape), nsets = length(names(xcape)), nintersects = NA)

