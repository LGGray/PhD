library(ggplot2)
library(ggrepel)

source('../../PhD/functions/edgeR.list.R')
load('../../datasets/XCI/chrX.Rdata')

edgeR <- edgeR.list('psuedobulk', filter=F)
names(edgeR) <- gsub('.edgeR-LRT', '', names(edgeR))

gene.set.enrichment <- function(list, gene.set, p.value){
    lapply(list, function(x){
        a <- length(intersect(x[x$FDR < p.value,1], gene.set))
        b <- length(intersect(x[x$FDR > p.value,1], gene.set))
        c <- length(setdiff(x[x$FDR < p.value,1], gene.set))
        d <- length(setdiff(x[x$FDR > p.value,1], gene.set))
        fisher.test(matrix(c(a,b,c,d), nrow=2), alternative='greater')
    })
}

result <- gene.set.enrichment(edgeR, rownames(chrX), p.value=0.05)

pvalue <- unlist(lapply(result, function(x) x$p.value))
odds.ratio <- unlist(lapply(result, function(x) x$estimate))
plot.data <- data.frame(celltype=names(result), pvalue=pvalue, odds.ratio=odds.ratio)

# dotplot x=odds.ratio, y=pvalue, color=pvalue, label=celltype
pdf('chrX.gene.set.enrichment.pdf')
ggplot(plot.data, aes(x=odds.ratio, y=-log10(pvalue), label=celltype)) +
    geom_point() +
    geom_text_repel() +
    scale_color_gradient(low='blue', high='red') +
    theme_bw() +
    theme(legend.position='none')
dev.off()