library(ggplot2)
library(reshape2)
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')

result <- read.delim('cpdb/analysis/significant_means.txt')

degs <- edgeR.list('psuedobulk', logfc=0.5)
names(degs) <- gsub('.edgeR-LRT', '', names(degs))

# subset result to only include names(degs)
lapply(strsplit(colnames(result)[-c(1:12)], '\\.'), function(x){
    names(degs) %in% x
})

lapply(names(degs), function(x){
    
}



colnames(result)
names(degs)