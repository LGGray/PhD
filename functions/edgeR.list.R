# Function to read in edgeR output files into list

deg.list <- function(path, method, logfc, filter=TRUE){
  if (method == NULL){
    stop("Please specify method")
  }
  files <- list.files(path, pattern=paste0(method,".txt"), full.names = T)
  deg.list <- lapply(files, function(x) read.delim(x, header=T, sep="\t"))
  names(deg.list) <- gsub('.txt', '', basename(files))
  if (filter != TRUE){
    return(deg.list)
    }
  if (filter == TRUE){
    deg.list <- lapply(deg.list, function(x) subset(x, FDR < 0.05 & abs(logFC) > logfc))
    deg.list <- deg.list[sapply(deg.list, function(x) nrow(x) > 0)]
    return(deg.list)
  }
}


