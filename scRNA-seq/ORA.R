# A function to perform over representation analysis (ORA) of MSigDB Hallmark gene sets
# in differentially expressed genes. This is done on upregulated and downregulated genes seperately
# and combined into a single output dataframe.

library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

ORA <- function(path, out, logFC=0){
  
  files <- list.files(path, 
                      pattern = ".txt")
  
  hallmark <- read.gmt('~/external/ClusterHome/gene.sets/h.all.v7.5.1.symbols.gmt')
  hallmark.list.up <- list()
  hallmark.list.down <- list()
  for ( file in files){
    infile <- read.delim(paste0(path, file))
    allOE_genes <- infile$gene
    FC.cuttoff <- logFC
    # Upregulated genes enrichment
    sigOE <- subset(infile, FDR < 0.05 & abs(logFC) > FC.cuttoff)
    sigOE_genes <- sigOE$gene
    
    up <- enricher(gene = sigOE_genes, 
                   universe = allOE_genes,
                   pAdjustMethod = "fdr", 
                   qvalueCutoff = 0.01,
                   TERM2GENE = hallmark)
    
    hallmark.list.up <- append(hallmark.list.up, up)
    names(hallmark.list.up)[length(hallmark.list.up)] <- gsub(".txt", '', file)
    
    # Downregulated genes enrichment
    sigOE <- subset(infile, FDR < 0.05 & abs(logFC) < FC.cuttoff)
    sigOE_genes <- sigOE$gene
    down <- enricher(gene = sigOE_genes, 
                     universe = allOE_genes,
                     pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.01,
                     TERM2GENE = hallmark)
    
    hallmark.list.down <- append(hallmark.list.down, down)
    names(hallmark.list.down)[length(hallmark.list.down)] <- gsub(".txt", '', file)
    
  }
  # Filter for empty list
  if(length(hallmark.list.up) > 0){
    hallmark.list.up <- hallmark.list.up[sapply(hallmark.list.up, function(x) length(x$ID) > 0)]
  }
  if(length(hallmark.list.down) > 0){
    hallmark.list.down <- hallmark.list.down[sapply(hallmark.list.down, function(x) length(x$ID) > 0)]
  }
  # Combine lists
  hallmark.list <- c(hallmark.list.up, hallmark.list.down)
  # Build list of dataframes and add logFC direction
  hallmark.df <- lapply(hallmark.list, data.frame)
  label <- c(rep('upregulated', length(hallmark.list.up)), 
             rep('downregulated', length(hallmark.list.down)))
  for(i in 1:length(hallmark.df)){hallmark.df[[i]]$`upregulated/downregulated` <- label[i]}
  # Extract data into single dataframe
  hallmark.file <- bind_rows(hallmark.df, .id = "column_label")
  # Write out to file
  write.table(hallmark.file, paste0(out, 'hallmark.txt') ,row.names=F, quote=F, sep="\t")
  return(list(hallmark.df, hallmark.file))
}

# Function to extract the genes and gene set from input gene list
search.ORA <- function(genes, infile){
  out <- list()
  for(i in 1:nrow(infile)){
    gene.list <- unique(unlist(strsplit(infile[i,9], "/")))
    xcape.list <- gene.list[gene.list %in% genes]
    if(length(xcape.list) > 0){
      out <- append(out, list(xcape.list))
      names(out)[length(out)] <- paste(infile[i,1], infile[i,2], sep='-')
    }
  }
  return(out)
}
