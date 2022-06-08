library(DOSE)
library(clusterProfiler)
#library(org.Hs.eg.db)
library(dplyr)

ORA <- function(path, out, logFC=0, pathway){
  
  files <- list.files(path, 
                      pattern = ".txt", full.names = T)
  # load("~/datasets/OneK1k/X_escape/escapees.Rdata")
  # load("~/datasets/OneK1k/X_escape/chrX.Rdata")
  
  gene.set <- read.gmt(pathway)
  gene.set.list.up <- list()
  gene.set.list.down <- list()
  for ( file in files){
    infile <- read.delim(file)
    allOE_genes <- infile$gene
    FC.cuttoff <- logFC
    # Upregulated genes enrichment
    sigOE <- subset(infile, FDR < 0.05 & abs(logFC) > FC.cuttoff)
    sigOE_genes <- sigOE$gene
    
    up <- enricher(gene = sigOE_genes, 
                   universe = allOE_genes,
                   pAdjustMethod = "fdr", 
                   qvalueCutoff = 0.01,
                   TERM2GENE = gene.set)
    
    gene.set.list.up <- append(gene.set.list.up, up)
    names(gene.set.list.up)[length(gene.set.list.up)] <- gsub(".txt", '', basename(file))
    
    # Downregulated genes enrichment
    sigOE <- subset(infile, FDR < 0.05 & abs(logFC) < FC.cuttoff)
    sigOE_genes <- sigOE$gene
    down <- enricher(gene = sigOE_genes, 
                     universe = allOE_genes,
                     pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.01,
                     TERM2GENE = gene.set)
    
    gene.set.list.down <- append(gene.set.list.down, down)
    names(gene.set.list.down)[length(gene.set.list.down)] <- gsub(".txt", '', basename(file))
    
  }
  # Filter for empty list
  if(length(gene.set.list.up) > 0){
    gene.set.list.up <- gene.set.list.up[sapply(gene.set.list.up, function(x) length(x$ID) > 0)]
  }
  if(length(gene.set.list.down) > 0){
    gene.set.list.down <- gene.set.list.down[sapply(gene.set.list.down, function(x) length(x$ID) > 0)]
  }
  # Combine lists
  gene.set.list <- c(gene.set.list.up, gene.set.list.down)
  # Build list of dataframes and add logFC direction
  gene.set.df <- lapply(gene.set.list, data.frame)
  label <- c(rep('upregulated', length(gene.set.list.up)), 
             rep('downregulated', length(gene.set.list.down)))
  for(i in 1:length(gene.set.df)){gene.set.df[[i]]$`upregulated/downregulated` <- label[i]}
  # Extract data into single dataframe
  gene.set.file <- bind_rows(gene.set.df, .id = "column_label")
  # Write out to file
  write.table(gene.set.file, out ,row.names=F, quote=F, sep="\t")
  return(gene.set.file)
}
# 
# xcape.gene.set.list <- list()
# index <- list()
# for(i in 1:nrow(result.file)){
#   gene.list <- unique(unlist(strsplit(result.file[i,9], "/")))
#   xcape.list <- gene.list[gene.list %in% rownames(escape)]
#   if(length(xcape.list) > 0){
#     xcape.gene.set.list <- append(xcape.gene.set.list, list(xcape.list))
#     names(xcape.gene.set.list)[length(xcape.gene.set.list)] <- result.file[i,2]
#   }
# }