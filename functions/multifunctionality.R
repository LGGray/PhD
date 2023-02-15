library(EGAD)

fibroblast <- read.delim('psuedobulk/Fibroblasts.edgeR-LRT.txt')
degs <- subset(fibroblast, abs(logFC) > 0.5 & FDR < 0.05)$gene
gene.set <- clusterProfiler::read.gmt('../../gene.sets/c5.go.bp.v7.5.1.symbols.gmt')
gene.set <- gene.set[,c(2,1)]

gene2annot <- subset(gene.set, gene %in% degs)
annotationlist <- unique(gene2annot$term)
annotations <- make_annotations(gene.set, degs, gene.set$term)
mf <- calculate_multifunc(annotations)
auc_mf <- auc_multifunc(annotations, mf$MF.rank)
hist <- plot_distribution(auc_mf, xlab="AUROC", med=FALSE, avg=FALSE)

mf <- mf[order(mf$MF.rank, decreasing=T),]
mf.chrX <- subset(mf, Gene %in% rownames(chrX))

