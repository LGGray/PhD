# function to return the most informative PC from PCA

calc.min.pc <- function(seurat.object){
  stdv <- seurat.object[['pca']]@stdev
  sum.stdv <- sum(seurat.object[['pca']]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  return(min.pc)
}