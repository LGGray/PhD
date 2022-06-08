library("UpSetR")

setwd('~/external/ClusterHome/datasets')
source('../R_code/functions/edgeR.list.R')
load('XCI/chrX.Rdata')
load('XCI/escapees.Rdata')

degs <- edgeR.list('GSE193770/edgeR-LRT/', logfc=0.5)

# Filter genes for X chromosome
Xchr <- lapply(degs, function(x){
  subset(x, gene %in% rownames(chrX))$gene
})
upset(fromList(Xchr), nsets = length(Xchr), nintersects = NA)
# Filter genes for XCI escape
xcape <- lapply(degs, function(x){
  subset(x, gene %in% rownames(escape))$gene
})
upset(fromList(xcape), nsets = length(xcape), nintersects = NA)

# Read in clusterProfiler files
files <- list.files(pattern='hallmark.txt', recursive = T,
                    full.names = T)
gene.set <- lapply(files, read.delim)
names(gene.set) <- c('UC', 'pSS', 'SLE', 'RA', 'MS')
gene.set.list <- lapply(gene.set, function(x){
 paste(x$column_label, x$ID, sep=':')
})

pdf('integrated/hallmark.upsetplot.pdf')
upset(fromList(gene.set.list), nsets=length(gene.set.list), nintersects = NA)
dev.off()

# Filter Terms for escape genes
gene.set.escape <- lapply(gene.set, function(x){
  x[gene.set.xcape(x),]})
gene.set.escape.list <- lapply(gene.set.escape, function(x){
  paste(x$column_label, x$ID, sep=':')
})
pdf('integrated/hallmark.escape.upsetplot.pdf')
upset(fromList(gene.set.escape.list), nsets=length(gene.set.escape.list), nintersects = NA)
dev.off()



