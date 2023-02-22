library(EGAD)

source('../../PhD/functions/edgeR.list.R')

# Read in GO gene sets
gene.set <- clusterProfiler::read.gmt('../../gene.sets/c5.go.v2022.1.Hs.symbols.gmt')[,c(2,1)]

# Create annotations
annotations <- make_annotations(gene.set, unique(gene.set$gene), unique(gene.set$term))
# Calculate multifunctionality
mf <- calculate_multifunc(annotations)
# Calculate AUC for Terms
auc_mf <- auc_multifunc(annotations, mf$MF.rank)
names(auc_mf) <- colnames(annotations)

# Read in DEGs
degs <- edgeR.list('psuedobulk', logfc=0.5)
# Match DEGs to MF score
degs.MF <- lapply(degs, function(x){
    mf[mf$Gene %in% x$gene,]
})

