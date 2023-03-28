library(gplots)
library(clusterProfiler)

source('../../PhD/functions/edgeR.list.R')
source('../../PhD/functions/chisq.test.degs.R')
source('../../PhD/functions/variance.list.R')
source('../../PhD/functions/chisq.test.var.R')

load('../../datasets/XCI/chrX.Rdata')

chrloc <- read.gmt('../../gene.sets/c1.all.v2023.1.Hs.symbols.gmt')
chrloc$chr <- paste0('chr', gsub("chr|p.+|q.+", "", chrloc$term))
chrloc <- subset(chrloc, !(chr %in% c('chrMT', 'chrY')))
chrloc.split <- split(chrloc, chrloc$chr)
chrloc.split$chrX <- data.frame(chr=rep('chrX', nrow(chrX)), gene=rownames(chrX))

# Calculate enrichment for DEGs
edgeR <- edgeR.list('psuedobulk', filter=FALSE)
result <- lapply(edgeR, function(x){
    res <- lapply(chrloc.split, function(y){
        chisq.test.degs(data=x, genes=y$gene, logfc=0.1)$p.value
    })
    names(res) <- names(chrloc.split)
    return(res)
})
result.df <- lapply(result, function(x){
    dplyr::bind_rows(x, .id='chr')
})

df <- dplyr::bind_rows(result.df, .id='cell')
df[is.na(df)] <- 1
pvals <- -log10(df[,-1])
rownames(pvals) <- df$cell

# Set a threshold for the p-value
pval_threshold <- 0.05
# Create a color palette
my_colors <- colorRampPalette(c("white", "red"))(100)
# Create a matrix of colors to use for the heatmap
colors <- matrix("white", nrow = nrow(pvals), ncol = ncol(pvals))
colors[pvals < pval_threshold] <- my_colors[100]

pdf('chromosome.enrichment.degs.pdf')
# Create a heatmap with color indicating significant p-value
heatmap.2(as.matrix(pvals), scale = "none", col = my_colors, trace = "none", key = TRUE,
          keysize = 1.5, density.info = "none", cexRow = 0.8, cexCol = 0.8,
          margins = c(10, 15), cellCol = colors)
dev.off()

# calculate enrichment for variance
var.list <- variance.list('variance', filter=F)
result <- lapply(var.list, function(x){
    res <- lapply(chrloc.split, function(y){
        chisq.test.var(data=x, genes=y$gene)$p.value
    })
    names(res) <- names(chrloc.split)
    return(res)
})
result.var <- lapply(result, function(x){
    dplyr::bind_rows(x, .id='chr')
})

# Plot heatmap of variance enrichment
df <- dplyr::bind_rows(result.var, .id='cell')
df[is.na(df)] <- 1
pvals <- -log10(df[,-1])
rownames(pvals) <- df$cell

# Set a threshold for the p-value
pval_threshold <- 0.05
# Create a color palette
my_colors <- colorRampPalette(c("white", "red"))(100)
# Create a matrix of colors to use for the heatmap
colors <- matrix("white", nrow = nrow(pvals), ncol = ncol(pvals))
colors[pvals < pval_threshold] <- my_colors[100]

pdf('chromosome.enrichment.var.pdf')
# Create a heatmap with color indicating significant p-value
heatmap.2(as.matrix(pvals), scale = "none", col = my_colors, trace = "none", key = TRUE,
          keysize = 1.5, density.info = "none", cexRow = 0.8, cexCol = 0.8,
          margins = c(10, 15), cellCol = colors)
dev.off()