# Function to perform fishers test of differentially expressed genes
# to determine enrichment in vector of genes
# The function expects edgeR output file

fisher.test.edgeR <- function(data, genes, logfc, direction='none')
{
  if (is.data.frame(data) == F) stop("data is not a data frame")
  if (is.vector(genes) == F) stop("genes is not a vector")
  if (is.numeric(logfc) == F) stop("logfc is not a number")

  if (direction == 'up') {
    a <- length(intersect(data[data$logFC > logfc & data$FDR < 0.05,'gene'], genes))
    b <- length(setdiff(data[data$logFC > logfc & data$FDR < 0.05,'gene'], genes))
    c <- length(intersect(data[data$FDR > 0.05,'gene'], genes))
    d <- length(setdiff(data[data$FDR > 0.05,'gene'], genes))
  } else if (direction == 'down') {
    a <- length(intersect(data[data$logFC < -logfc & data$FDR < 0.05,'gene'], genes))
    b <- length(setdiff(data[data$logFC < -logfc & data$FDR < 0.05,'gene'], genes))
    c <- length(intersect(data[data$FDR > 0.05,'gene'], genes))
    d <- length(setdiff(data[data$FDR > 0.05,'gene'], genes))
  } else {
    a <- length(intersect(data[abs(data$logFC) > logfc & data$FDR < 0.05,'gene'], genes))
    b <- length(setdiff(data[abs(data$logFC) > logfc & data$FDR < 0.05,'gene'], genes))
    c <- length(intersect(data[data$FDR > 0.05,'gene'], genes))
    d <- length(setdiff(data[data$FDR > 0.05,'gene'], genes))
  }
  
  test = fisher.test(matrix(c(a,b,c,d), nrow=2), alternative='greater')
  return(test)
}