# Function to match TF targetting eGenes to the location of eQTL in TFBS
# targets.file = The gene targets of a set TFs
# TFBS.intersect = Output of a list of eQTLs overlapping remap CRM TFBS. Output of bedtools intersect

library(dplyr)
TFBS.eQTL <- function(targets.file, TFBS.intersect){
  result.list <- list()
  for(i in 1:nrow(targets.file)){
    tmp <- subset(TFBS.intersect, V14 %in% targets.file[i,1] & V13 %in% targets.file[i,3])
    if(nrow(tmp) > 0){
      TF.list <- lapply(tmp$V4, function(x) unlist(strsplit(x, ',')))
      foo <- lapply(TF.list, function(x) grep(paste0(targets.file[i,2], '\\b'), x))
      if(length(foo[[1]]) > 0){
        bar <- cbind(tmp[,c(1,2,3,10,11,12)], targets.file[i,1:3])
        colnames(bar) <- c('TFBS.chr', 'TFBS.start', 'TFBS.end', 'eQTL.chr', 
                           'eQTL.start', 'eQTL.end', 'cell_type', 'TF', 'eGene')
        result.list[[length(result.list)+1]] <- bar
      }
    }
  }
  result.df <- bind_rows(result.list, .id=NULL)
  return(result.df)
}


