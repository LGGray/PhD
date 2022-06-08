
library(dplyr)
library(org.Hs.eg.db)

magma.infile <- function(deg.list, direction=NULL, outfile){
  # if (direction != 'up'| direction != 'down') stop("direction must be 'up' or 'down'")
  if(direction == 'up'){
    degs <- lapply(deg.list, function(x) subset(x, logFC > 0)[,1])
  }
  if(direction == 'down'){
    degs <- lapply(deg.list, function(x) subset(x, logFC < 0)[,1])
  }
  names(degs) <- names(deg.list)
  degs.df <- lapply(degs, data.frame)
  degs.file <- bind_rows(degs.df, .id = "column_label")
  colnames(degs.file) <- c('celltype', 'SYMBOL')
  
  converted <- AnnotationDbi::select(org.Hs.eg.db, degs.file$SYMBOL, 
                                     columns='ENTREZID', keytype = 'SYMBOL')
  converted <- converted[!is.na(converted$ENTREZID),]
  degs.file <- merge(converted, degs.file, by='SYMBOL')[,c(3,1,2)]
  degs.file <- degs.file[sort.int(degs.file$celltype, index.return=T)$ix,]
  
  write.table(degs.file, outfile, 
              row.names=F, quote=F, sep="\t")
}

magma.infile(SDY998, direction = 'down', '~/external/ClusterHome/datasets/SDY998/magma/downreg.txt')




