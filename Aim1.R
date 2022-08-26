# Code to produce figures for PhD Aim 1

setwd('~/external/ClusterHome/datasets/integrated/')

library(ggplot2)
library(reshape2)
library(UpSetR)
library(clusterProfiler)
library(combinat)
library(gplots)
source('../../R_code/functions/edgeR.list.R')
source('../../R_code/functions/chisq.test.degs.R')
load('../XCI/escapees.Rdata')
load('../XCI/chrX.Rdata')
source('../../R_code/colour.dictionary.R')

# All genes
all <- edgeR.list('psuedobulk/', filter=F)
all <- all[-17]
# Differentially expressed genes
deg <- edgeR.list('psuedobulk/', logfc=0.5)
deg <- deg[-17]
# Those downregulated i.e more abundant in disease
down <- lapply(deg, function(x) subset(x, logFC < 0))
# Those upregulated i.e more abundant in controls
up <- lapply(deg, function(x) subset(x, logFC > 0))

# Average DEGs across celltypes and within up and down regulated sets
mean(unlist(lapply(deg, nrow)))
mean(unlist(lapply(up, nrow)))
mean(unlist(lapply(down, nrow)))

# Count number of total and escape DEGs in up/downregulated sets
deg.count <- lapply(deg, function(x){
  up <- subset(x, logFC > 0)
  down <- subset(x, logFC < 0)
  df.up <- data.frame(escape=length(up$gene[up$gene %in% rownames(chrX)]),total=length(up$gene), direction='Up')
  df.down <- data.frame(escape=length(down$gene[down$gene %in% rownames(chrX)]), total=length(down$gene), direction='Down')
  df <- rbind(df.up, df.down)
  df$total <- df$total-df$escape
  return(df)
})
deg.count <- dplyr::bind_rows(deg.count, .id='celltype')
# Barplot
ggplot(melt(deg.count), aes(x=celltype, y=value, fill=direction)) + 
  geom_bar(stat='identity', position = position_dodge()) + 
  scale_fill_brewer(palette='Paired') +
  theme(axis.text.y = element_text(size=14)) +
  coord_flip() +
  xlab('') + ylab('DEG count') + labs(fill='Direction')

# Calculate enrichment of chrX genes
enrichment <- lapply(all, function(x){
  up <- subset(x, logFC > 0)
  down <- subset(x, logFC < 0)
  data.frame(up=chisq.test.degs(up, rownames(chrX), logfc=0.5)$p.value, 
             down=chisq.test.degs(down, rownames(chrX), logfc=0.5)$p.value)
})
enrichment <- dplyr::bind_rows(enrichment, .id='celltype')
enrichment <- enrichment[enrichment$down < 0.05,]

# Overlap of chrX between celltypes with chi enrichment
cells <- enrichment$celltype
lst <- lapply(down[cells], function(x){
  subset(x, gene %in% rownames(chrX))$gene
})
upset.data <- fromList(lst)
rownames(upset.data) <- unique(unlist(lst))
pdf('chrX.upsetplot.pdf')
upset(upset.data, nsets = ncol(upset.data), nintersects = NA,
      sets=colnames(upset.data), keep.order=T, sets.bar.color= cell.colourdict[colnames(upset.data)], order.by = 'freq')
dev.off()
upset.data[order(rowSums(upset.data), decreasing = T),]

# Load Gene Set
gene.set <- read.gmt('../../gene.sets/c3.tft.gtrd.v7.5.1.symbols.gmt')

# Over-representation analysis for downregulated genes in enriched cells
ora <- lapply(all, function(x){
  enricher(gene = subset(x, FDR < 0.05 & logFC < -0.5)$gene, 
           universe = subset(x, logFC < -0.5)$gene,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.01,
           TERM2GENE = gene.set)
})
# Remove empty list
ora <- ora[unlist(lapply(ora, function(x) !is.null(x)))]

# Ceate heatmap for all genes
ora.result <- lapply(ora, function(x) x@result[,c(1,7,9)])
ora.result <- dplyr::bind_rows(ora.result, .id='celltype')
# Filter
ora.result <- ora.result[ora.result$qvalue < 0.05,]
# Clustered heatmap
# m <- tidyr::pivot_wider(ora.result, names_from = 'celltype', values_from = 'qvalue')
# m[is.na(m)] <- 1
# m <- as.matrix(m[,-1])
# # m <- m[rowSums(m) < 7,]
# # m <- m[,colSums(m) < nrow(m)]
# 
# clust <- hclust(dist(t(m)))
# ggplot(ora.result, aes(x=celltype, y=ID, fill=qvalue)) +
#   geom_tile() +
#   scale_fill_gradient(low = "red", high = "white", limits=c(0,1), breaks = c(0, 0.5, 1)) +
#   scale_x_discrete(limits = colnames(m)[clust$order]) +
#   theme(axis.text.x = element_text(angle=45, hjust=1))

pdf('clusterProfiler/chrX/GTRD.down.pdf', height = 10)
ggplot(ora.result, aes(x=celltype, y=ID, colour=qvalue, size=Count)) + 
  geom_point() +
  scale_colour_gradient(low = "red", high = "purple") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  xlab('') + ylab('')
dev.off()

unique(ora.result$ID)[unique(ora.result$ID) %in% rownames(chrX)]
# Subset for terms containing escape genes
ora.x <- lapply(ora, function(x) x@result[,c(1,7,9,8)])
ora.x <- dplyr::bind_rows(ora.x, .id='celltype')
# Filter
ora.x <- subset(ora.x, qvalue < 0.05)
lst <- strsplit(ora.x$geneID, '/')
lst <- lapply(lst, function(x) x[x %in% rownames(chrX)])
ora.x <- ora.x[sapply(lst, function(x) length(x) > 0),1:4]



pdf('clusterProfiler/chrX/GTRD.down.chrX.pdf',height = 10)
ggplot(ora.x, aes(x=celltype, y=ID, colour=qvalue, size=Count)) + 
  geom_point() +
  scale_colour_gradient(low = "red", high = "purple") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  xlab('') + ylab('')
dev.off()


# # Clustered heatmap
# m <- tidyr::pivot_wider(ora.x, names_from = 'celltype', values_from = 'qvalue')
# m <- as.matrix(m[,-1])
# m[is.na(m)] <- 1
# clust <- hclust(dist(t(m)))
# ggplot(ora.x, aes(x=celltype, y=ID, fill=qvalue)) +
#   geom_tile() +
#   scale_fill_gradient(low = "red", high = "white", limits=c(0,1), breaks = c(0, 0.5, 1)) +
#   scale_x_discrete(limits = colnames(m)[clust$order]) +
#   theme(axis.text.x = element_text(angle=45, hjust=1))

# Over-representation analysis for downregulated genes
ora.down <- lapply(down, function(x){
  enricher(gene = subset(x, FDR < 0.05 & abs(logFC) > 0.5)$gene, 
           universe = x$gene,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.01,
           TERM2GENE = gene.set)
})
# Remove empty list
ora.down <- ora.down[unlist(lapply(ora.down, function(x) !is.null(x)))]

# Ceate heatmap for all genes
ora.result <- lapply(ora.down, function(x) x@result[,c(1,7)])
ora.result <- dplyr::bind_rows(ora.result, .id='celltype')
ora.result <- subset(ora.result, qvalue < 0.05)
# Clustered heatmap
m <- tidyr::pivot_wider(ora.result, names_from = 'celltype', values_from = 'qvalue')
m <- as.matrix(m[,-1])
clust <- hclust(dist(t(m)))
ggplot(ora.result, aes(x=celltype, y=ID, fill=qvalue)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "white") +
  scale_x_discrete(limits = colnames(m)[clust$order]) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

# Subset for terms containing escape genes
ora.xcape <- lapply(ora.down, function(x) x@result[,c(1,7,8)])
ora.xcape <- dplyr::bind_rows(ora.xcape, .id='celltype')
ora.xcape <- subset(ora.xcape, qvalue < 0.05)
lst <- strsplit(ora.xcape$geneID, '/')
lst <- lapply(lst, function(x) x[x %in% rownames(escape)])
ora.xcape <- ora.xcape[sapply(lst, function(x) length(x) > 0),1:3]
# Clustered heatmap
m <- tidyr::pivot_wider(ora.xcape, names_from = 'celltype', values_from = 'qvalue')
m <- as.matrix(m[,-1])
m[is.na(m)] <- 1
clust <- hclust(dist(t(m)))
ggplot(ora.xcape, aes(x=celltype, y=ID, fill=qvalue)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "white") +
  scale_x_discrete(limits = colnames(m)[clust$order]) +
  theme(axis.text.x = element_text(angle=45, hjust=1))


# correlation <- lapply(combn(1:length(deg), 2, simplify = F), function(x){
#   df <- merge(deg[[x[1]]][,1:2], deg[[x[2]]][,1:2], by='gene')
#   data.frame(combo=paste0(names(deg)[x[1]],':',names(deg)[x[2]]), r=cor(df$logFC.x, df$logFC.y, method='pearson'))
# })
# correlation <- dplyr::bind_rows(correlation, .id=NULL)
# 
# subset(correlation, abs(r) > 0.8)
# 
# test <- merge(deg[['B_intermediate']][,1:2], deg[['B_memory']][,1:2], by='gene')
# cor(test$logFC.x, test$logFC.y, method='pearson')
# plot(test$logFC.x, test$logFC.y)
# abline(lm(test$logFC.x ~ test$logFC.y), col='red')


# Upset plot
xcape <- lapply(deg, function(x){
  x$gene[x$gene %in% rownames(escape)]
})
xcape.overlap <- fromList(xcape)
rownames(xcape.overlap) <- unique(unlist(xcape))
upset(xcape.overlap, nsets=length(xcape))
xcape.overlap$count <- rowSums(xcape.overlap)
xcape.overlap <- xcape.overlap[order(xcape.overlap$count, decreasing=T),]

# Functional analysis of escape genes
immport <- read.delim('../../gene.sets/Immport.GeneList.txt')
xcape.import <- lapply(xcape, function(x){
  subset(immport, Symbol %in% x)
})
xcape.import

immune <- read.delim('../XCI/X.immune.txt')
xcape.immune <- lapply(xcape, function(x){
  subset(immune, gene %in% x)
})
xcape.immune <- xcape.immune[sapply(xcape.immune, function(x) nrow(x) > 0)]
xcape.immune <- dplyr::bind_rows(xcape.immune, .id='celltype')
write.table(xcape.immune, 'xcape.immune.txt', sep='\t', row.names=F, quote=F)

gene.set <- read.gmt('../../gene.sets/h.all.v7.5.1.symbols.gmt')
xcape.geneset <- lapply(xcape[enriched], function(x){
  subset(gene.set, gene %in% x)
})




