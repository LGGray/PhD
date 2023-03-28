# Chi-squared test of enrichmnet for variance output files
chisq.test.var <- function(data, genes)
{
  if (is.data.frame(data) == F) stop("data is not a data frame")
  if (is.vector(genes) == F) stop("genes is not a vector")
  a <- length(intersect(data[data$FDR < 0.05,'gene'], genes))
  b <- length(setdiff(data[data$FDR < 0.05,'gene'], genes))
  c <- length(intersect(data[data$FDR > 0.05,'gene'], genes))
  d <- length(setdiff(data[data$FDR > 0.05,'gene'], genes))
  test = (chisq.test(matrix(c(a,b,c,d), nrow=2)))
  return(test)
}