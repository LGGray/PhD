binary.matrix <- function(deg.list, mode='both'){
  if(mode=='both'){
    genes <- unique(unlist(lapply(deg.list, function(x){
      x$gene
    })))
    mtx <- matrix(nrow = length(genes), ncol=length(deg.list))
    rownames(mtx) <- genes
    colnames(mtx) <- names(deg.list)                    
    for( cell in colnames(mtx)){
      mtx[,cell] <- ifelse(rownames(mtx) %in% deg.list[[cell]]$gene, 1, 0)
    }
    return(mtx)
  } 
  if(mode == 'up'){
    genes <- unique(unlist(lapply(deg.list, function(x){
      subset(x, logFC > 0)$gene
    })))
    mtx <- matrix(nrow = length(genes), ncol=length(deg.list))
    rownames(mtx) <- genes
    colnames(mtx) <- names(deg.list)                    
    for( cell in colnames(mtx)){
      mtx[,cell] <- ifelse(rownames(mtx) %in% deg.list[[cell]]$gene, 1, 0)
    }
    return(mtx)
  }
  if(mode == 'down'){
    genes <- unique(unlist(lapply(deg.list, function(x){
      subset(x, logFC < 0)$gene
    })))
    mtx <- matrix(nrow = length(genes), ncol=length(deg.list))
    rownames(mtx) <- genes
    colnames(mtx) <- names(deg.list)                    
    for( cell in colnames(mtx)){
      mtx[,cell] <- ifelse(rownames(mtx) %in% deg.list[[cell]]$gene, 1, 0)
    }
    return(mtx)
  }
}


