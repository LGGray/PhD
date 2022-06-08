# Function to perform Chi-squared test of differentially expressed genes
# to determine enrichment in vector of genes
# The function expects edgeR output file

chisq.test.degs <- function(infile, genes, logFC)
{
  if (is.vector(genes) == F) stop("genes is not a vector")
  if (is.numeric(logFC) == F) stop("logFC is not a number")
  total.genes <- length(infile$gene)
  total.list <- length(which(infile$gene %in% genes))
  deg.total <- length(which(infile$FDR < 0.05 & abs(infile$logFC) > logFC))
  deg.list <- length(which(infile$FDR < 0.05 & abs(infile$logFC) > logFC & infile$gene %in% genes))
  ctable <- matrix(c(deg.list, deg.total-deg.list, total.list-deg.list, 
                     (total.genes-total.list)-(deg.total-deg.list)), nrow=2)
  test = (chisq.test(ctable))
  return(test)
}