# Analysing the overlap between the diseases
library(tidygraph)
library(tidyverse)
library(ggraph)
library(ggvenn)
setwd('~/external/ClusterHome/datasets/')
load('OneK1k/X_escape/escapees.Rdata')

# Read in files
files <- list.files(pattern='kegg.txt', full.names=T, recursive = T)
gene.set <- lapply(files, read.delim)
names(gene.set) <- c('UC', 'pSS', 'SLE', 'RA')
gene.set <- lapply(gene.set, function(x) x[!is.na(x[,8]),])
gene.set <- lapply(gene.set, function(x) subset(x, x[,8] < 0.01))
# Add new column of cell:term
gene.set <- lapply(gene.set, function(x) x <- cbind(x, 
                                                    term=paste(x[,1],x[,3], sep=':')))

# venn diagram
terms <- list('UC'=gene.set[[1]]$term, 'pSS'=gene.set[[2]]$term,
     'SLE'=gene.set[[3]]$term, 'RA'=gene.set[[4]]$term)
ggvenn(terms)

# Export cytoscape file for all gene sets
cytoscape <- lapply(gene.set, function(x) x[,c(12,8)])
cytoscape <- dplyr::bind_rows(cytoscape, .id = 'column_name')
write.table(cytoscape, 'integrated/figures/kegg.cytoscape.txt', row.names=F, quote=F, sep="\t")

# Filter for gene sets containing escapers and produce cytoscape file
gene.set.xcape <- function(df){
  index <- list()
  xcape.gene.set.list <- list()
  for(i in 1:nrow(df)){
    gene.list <- unique(unlist(strsplit(df[i,9], "/")))
    xcape.list <- gene.list[gene.list %in% rownames(escape)]
    if(length(xcape.list) > 0){
      index <- append(index, i)
      xcape.gene.set.list <- append(xcape.gene.set.list, list(xcape.list))
      names(xcape.gene.set.list)[length(xcape.gene.set.list)] <- df[i,12]
    }
  }
  return(unlist(index))
}

escape.cytoscape <- lapply(gene.set, function(x) x[gene.set.xcape(gene.set[[1]]),c(12,8)])
escape.cytoscape <- lapply(escape.cytoscape, function(x) x[!is.na(x$term),])
escape.cytoscape.wide <- dplyr::bind_rows(escape.cytoscape, .id = 'column_name')
write.table(escape.cytoscape.wide, 'integrated/figures/kegg.escape.cytoscape.txt', row.names=F, quote=F, sep="\t")

# Calculate pearson correlation of matching terms
ORA.cor <- function(x, y){
  common <- intersect(x$term, y$term)
  if(length(common)==0){
    return('No common features')
  }else{
    a <- subset(x, term %in% common)
    z <- subset(y, term %in% common)
    print(dim(a))
    return(cor(a[,2], z[,2], method='pearson'))
  }
}

test.cor <- function(file.list){
  # Vector of combinations from 1:4
  comb <- data.frame(gtools::combinations(length(file.list), 2, 1:length(file.list)))
  for(i in 1:nrow(comb)){
    r <- ORA.cor(file.list[[comb[i,1]]], file.list[[comb[i,2]]])
    if(is.numeric(r) & !is.na(r) & r!="No common features"){
      print(paste(names(gene.set)[comb[i,1]], names(gene.set)[comb[i,2]],
                  r, sep = ':'))
    }
  }
}

test.cor(escape.cytoscape)


