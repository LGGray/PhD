library("UpSetR")

setwd('~/external/ClusterHome/datasets')
source('../R_code/functions/edgeR.list.R')
source('~/PhD/colour.dictionary.R')
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

# Compare enrichment across disease and cell type
filename <- 'reactome'
files <- list.files(pattern=paste0(filename, '.txt'), recursive = T,
                    full.names = T)
gene.set <- lapply(files, read.delim)
names(gene.set) <- c('UC', 'pSS', 'MS', 'SLE', 'RA')
gene.set.list <- lapply(gene.set, function(x){
 paste(x$column_label, x$ID, sep=':')
})
# Build upsetplot
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

# Compare enrichment across cell types within disease
study <- 'GSE193770' # Change disease below
set <- 'kegg'

gene.set <- read.delim(paste0(study,'/','clusterProfiler','/',set,'.txt'))
gene.set$column_label <- factor(gene.set$column_label)
gene.set.list <- split(gene.set, gene.set$column_label)
gene.set.list <- lapply(gene.set.list, function(x) x$ID)

pdf(paste0(study, '/', set, '.upsetplot.pdf'))
data <- fromList(gene.set.list)
upset(data, nsets = ncol(data), nintersects = NA,
      sets=colnames(data), keep.order=T, 
      sets.bar.color=cell.colourdict[colnames(data)],
      main.bar.color = study.colourdict['MS'])
dev.off()

# Compare enrichment across disease in the same cell type
filename <- 'reactome'
files <- list.files(pattern=paste0(filename, '.txt'), recursive = T,
                    full.names = T)
gene.set <- lapply(files, read.delim)
names(gene.set) <- c('UC', 'pSS', 'MS', 'SLE', 'RA')
# Find cells represented across data
cells <- unique(unlist(lapply(gene.set, function(x) x$column_label)))
dir <- ('integrated/upsetplots/cells')
for(i in 1:length(cells)){
  gene.set.list <- lapply(gene.set, function(x){
    subset(x, column_label %in% cells[i])$ID
  })
  data <- fromList(gene.set.list)
  pdf(paste0(dir, '/', filename, '.', cells[i], '.pdf'))
  upset(data, nsets = ncol(data), nintersects = NA,
        sets=colnames(data), keep.order=T,
        sets.bar.color=unlist(study.colourdict[colnames(data)]),
        main.bar.color = cell.colourdict[cells[i]])
  dev.off()
}

cell=cells[23]
gene.set.list <- lapply(gene.set, function(x){
  subset(x, column_label %in% cell)$ID
})
data <- fromList(gene.set.list)
pdf(paste0(dir, '/', filename, '.', cell, '.pdf'))
upset(data, nsets = ncol(data), nintersects = NA,
      sets=colnames(data), keep.order=T,
      sets.bar.color=unlist(study.colourdict[colnames(data)]),
      main.bar.color = cell.colourdict[cell])
dev.off()

