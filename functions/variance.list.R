# function to read in variance files into a list

variance.list <- function(path, filter=TRUE){
  files <- list.files(path, pattern=".txt", full.names = T)
  var.list <- lapply(files, function(x) read.delim(x, header=T, sep="\t"))
  names(var.list) <- gsub('.txt', '', basename(files))
  if (filter != TRUE){
    return(var.list)
    }
  if (filter == TRUE){
    var.list <- lapply(var.list, function(x) subset(x, FDR < 0.05))
    var.list <- var.list[sapply(var.list, function(x) nrow(x) > 0)]
    return(var.list)
  }
}