library(chromoMap)
library(biomaRt)
library(ggplot2)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") # for GRCh37 version = 'GRCh37'

source('~/external/ClusterHome/R_code/functions/edgeR.list.R')
load('~/external/ClusterHome/datasets/XCI/chrX.Rdata')

setwd('~/external/ClusterHome/datasets/integrated/')

up <- edgeR.list('~/external/ClusterHome/datasets/integrated/psuedobulk/upregulated/', logfc=0.5)
up.X <- lapply(up, function(x){
  subset(x, gene %in% rownames(chrX))$gene
})
up.lst <- unique(unlist(up.X))

down <- edgeR.list('~/external/ClusterHome/datasets/integrated/psuedobulk/downregulated/', logfc=0.5)
down.X <- lapply(down, function(x){
  subset(x, gene %in% rownames(chrX))$gene
})
down.lst <- unique(unlist(down.X))

up.coords <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
      filters = 'external_gene_name',
      values = up.lst,
      mart = ensembl)

up.coords$arm <- ifelse(up.coords$start_position < 60509062 & up.coords$end_position < 60509062, 'p', 'q')
up.split <- split(up.coords, up.coords$arm)

down.coords <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
                   filters = 'external_gene_name',
                   values = down.lst,
                   mart = ensembl)
down.coords$arm <- ifelse(down.coords$start_position < 60509062 & down.coords$end_position < 60509062, 'p', 'q')
down.split <- split(down.coords, down.coords$arm)

wilcox.test(sort(up.split$p$start_position), sort(down.split$p$start_position), alternative = 'less')
wilcox.test(up.split$q$start_position, down.split$q$start_position)

plot.data <- data.frame(up=up.split$p$start_position, down=c(down.split$p$start_position, rep(NA, 11)))
ggplot(melt(plot.data), aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9)

boxplot(up.split$p$start_position, down.split$p$start_position)
boxplot(up.split$q$start_position, down.split$q$start_position)

chr_file <- data.frame(V1='X', V2=10001, V3=156030895, V4=60509062)
names(chr_file) <- NULL

write.table(coords, 'annotation.txt', quote=F,
            sep="\t",row.names=FALSE,col.names=FALSE)
write.table(down.coords, 'down.annotation.txt', quote=F,
            sep="\t",row.names=FALSE,col.names=FALSE)
write.table(chr_file, 'chr_file.txt', quote=F,
            sep="\t",row.names=FALSE,col.names=FALSE)


chromoMap('chr_file.txt', 'annotation.txt')
