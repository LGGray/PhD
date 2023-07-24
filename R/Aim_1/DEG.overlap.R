
source('../PhD/functions/edgeR.list.R')
load('../datasets/XCI/chrX.Rdata')

pSS <- deg.list('pSS_GSE157278/differential.expression/edgeR/', logfc=0.5)
pSS <- lapply(pSS, function(x) x[x$gene %in% rownames(chrX),])
MS <- deg.list('MS_GSE193770/differential.expression/edgeR/', logfc=0.5)
MS <- lapply(MS, function(x) x[x$gene %in% rownames(chrX),])
UC <- deg.list('UC_GSE125527/differential.expression/edgeR/', logfc=0.5)
UC <- lapply(UC, function(x) x[x$gene %in% rownames(chrX),])
CD_colon <- deg.list('CD_Kong/colon/differential.expression/edgeR/', logfc=0.5)
CD_colon <- lapply(CD_colon, function(x) x[x$gene %in% rownames(chrX),])
CD_TI <- deg.list('CD_Kong/TI/differential.expression/edgeR/', logfc=0.5)
CD_TI <- lapply(CD_TI, function(x) x[x$gene %in% rownames(chrX),])
SLE <- deg.list('lupus_Chun/differential.expression/edgeR/', logfc=0.5)
SLE <- lapply(SLE, function(x) x[x$gene %in% rownames(chrX),])

# Find common celltypes as names of the lists
common <- Reduce(intersect, list(names(pSS), names(UC), names(CD_colon), names(CD_TI), names(SLE)))

pSS <- pSS[common]
UC <- UC[common]
CD_colon <- CD_colon[common]
CD_TI <- CD_TI[common]
SLE <- SLE[common]

# Create heatmap of chrX expression in each celltype and disease
common_genes <- lapply(common, function(celltype) {
  unique(unlist(list(pSS[[celltype]]$gene, UC[[celltype]]$gene, CD_colon[[celltype]]$gene, CD_TI[[celltype]]$gene, SLE[[celltype]]$gene)))
})
names(common_genes) <- common

matrices <- lapply(common, function(celltype) {
    genes <- common_genes[[celltype]]
    mtx <- matrix(0, nrow=length(common_genes[[celltype]]), ncol=5)
    rownames(mtx) <- genes
    colnames(mtx) <- c('pSS', 'UC', 'CD_colon', 'CD_TI', 'SLE')
    mtx[,1] <- pSS[[celltype]]$logFC[match(genes, pSS[[celltype]]$gene)]
    mtx[,2] <- UC[[celltype]]$logFC[match(genes, UC[[celltype]]$gene)]
    mtx[,3] <- CD_colon[[celltype]]$logFC[match(genes, CD_colon[[celltype]]$gene)]
    mtx[,4] <- CD_TI[[celltype]]$logFC[match(genes, CD_TI[[celltype]]$gene)]
    mtx[,5] <- SLE[[celltype]]$logFC[match(genes, SLE[[celltype]]$gene)]
    return(mtx)
})

pdf('Tem_Trm_cytotoxic_T_cells.chrX.heatmap.pdf', width=10, height=11)
ggplot(melt(matrices[[7]]), aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low='blue', high='red') + 
  theme_bw() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  ylab("Gene") + xlab("Cell type") + ggtitle("Tem/Trm cytotoxic T cells: chrX logFC ") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

library(gplots)

# Define a function to create a heatmap for a given matrix
create_heatmap <- function(mtx) {
  heatmap.2(mtx, trace='none', col=rev(colorRampPalette(c('blue', 'white', 'red'))(100)), Rowv=FALSE, Colv=FALSE, margins=c(10,10))
}

pdf('DC2.chrX.heatmap.pdf')
heatmap.2(matrices[[1]], trace='none', col=rev(colorRampPalette(c("blue", "white", "red"))(100)), 
          scale='none', margins=c(10,10), main='DC2 chrX genes', xlab='', ylab='',
          srtCol=45, key=FALSE, dendrogram='none')
dev.off()

# Create a heatmap for each cell type
for (celltype in common) {
  mtx <- matrices[[celltype]]
  create_heatmap(mtx)
}