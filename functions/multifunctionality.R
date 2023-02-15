library(EGAD)

fibroblast <- read.delim('psuedobulk/Fibroblasts.edgeR-LRT.txt')
degs <- subset(fibroblast, abs(logFC) > 0.5 & FDR < 0.05)$gene
gene.set <- clusterProfiler::read.gmt('../../gene.sets/c5.go.v2022.1.Hs.symbols.gmt')[,c(2,1)]

annotations <- make_annotations(gene.set, unique(gene.set$gene), unique(gene.set$term))
mf <- calculate_multifunc(annotations)
auc_mf <- auc_multifunc(annotations, mf$MF.rank)
pdf('AUC_mf.pdf')
hist <- plot_distribution(auc_mf, xlab="AUROC", med=FALSE, avg=FALSE)
dev.off() 

subset(mf, Gene %in% degs)