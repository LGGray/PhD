# Function to perform fishers test of differentially expressed genes
# to determine enrichment in vector of genes
# The function expects edgeR output file

fishers.test.degs <- function(infile, genes, logfc)
{
  if (is.vector(genes) == F) stop("genes is not a vector")
  if (is.numeric(logfc) == F) stop("logfc is not a number")
  total.genes <- length(infile$gene)
  total.list <- length(which(infile$gene %in% genes))
  a <- length(infile[subset(infile, FDR < 0.05 & abs(logFC) > logfc)$gene %in% genes,1])
  b <- length(infile[!subset(infile, FDR < 0.05 & abs(logFC) > logfc)$gene %in% genes,1])
  c <- length(infile[subset(infile, FDR > 0.05)$gene %in% genes, 1])
  d <- length(infile[!subset(infile, FDR > 0.05)$gene %in% genes, 1])
  ctable <- matrix(c(a,b,c,d), nrow=2)
  test = fisher.test(ctable)
  return(test)
}


