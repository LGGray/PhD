# Function to read in edgeR output files into list

deg.list <- function(path, logfc, filter=TRUE){
  files <- list.files(path, full.names = T)
  deg.list <- lapply(files, function(x) read.delim(x, header=T, sep="\t"))
  names(deg.list) <- gsub('.txt', '', basename(files))
  if (filter != TRUE){
    return(deg.list)
    }
  if (filter == TRUE){
    deg.list <- lapply(deg.list, function(x) subset(x, FDR < 0.05 & abs(logFC) > logfc))
    return(deg.list)
  }
}


