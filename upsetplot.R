library("UpSetR")

setwd('~/external/ClusterHome/datasets')
source('../R_code/functions/edgeR.list.R')
load('XCI/chrX.Rdata')
load('XCI/escapees.Rdata')

directory <- 'SDY998'
degs <- edgeR.list(paste0(directory,'/edgeR-LRT/'), logfc=0.5)

# Filter genes for X chromosome
Xchr <- lapply(degs, function(x){
  subset(x, gene %in% rownames(chrX))$gene
})
pdf(paste0(directory, '/chrX.upsetplot.pdf'))
data <- fromList(Xchr)
upset(data, nsets = ncol(data), nintersects = NA,
      sets=colnames(data), keep.order=T, sets.bar.color= cell.colourdict[colnames(data)])
dev.off()

# Filter genes for XCI escape
xcape <- lapply(degs, function(x){
  subset(x, gene %in% rownames(escape))$gene
})
pdf(paste0(directory, '/escape.upsetplot.pdf'))
data <- fromList(xcape)
upset(data, nsets = ncol(data), nintersects = NA,
      sets=colnames(data), keep.order=T, sets.bar.color= cell.colourdict[colnames(data)])
dev.off()

# Read in clusterProfiler files
filename <- 'reactome'
files <- list.files(pattern=paste0(filename, '.txt'), recursive = T,
                    full.names = T)
gene.set <- lapply(files, read.delim)
names(gene.set) <- c('UC', 'pSS', 'SLE', 'RA', 'MS')
gene.set.list <- lapply(gene.set, function(x){
 paste(x$column_label, x$ID, sep=':')
})

pdf(paste0('integrated/', filename, '.upsetplot.pdf'))
data <- fromList(gene.set.list)
upset(data, nsets = ncol(data), nintersects = NA,
      sets=colnames(data), keep.order=T, sets.bar.color=study.colourdict[colnames(data)])
dev.off()

# Filter Terms for escape genes
gene.set.escape <- lapply(gene.set, function(x){
  x[gene.set.xcape(x),]})
gene.set.escape.list <- lapply(gene.set.escape, function(x){
  paste(x$column_label, x$ID, sep=':')
})
pdf(paste0('integrated/', filename, '.escape.upsetplot.pdf'))
data <- fromList(gene.set.escape.list)
upset(data, nsets = ncol(data), nintersects = NA,
      sets=colnames(data), keep.order=T, sets.bar.color=study.colourdict[colnames(data)])
dev.off()



