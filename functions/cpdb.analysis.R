library(ggplot2)
library(reshape2)
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
load('../../datasets/XCI/chrX.Rdata')

interactions <- read.delim('cpdb/degs_analysis/relevant_interactions.txt')
interactions[,1:11]

unlist(strsplit(interactions$interacting_pair, '_')) %in% rownames(chrX)

cell.interactions <- lapply(1:nrow(interactions), function(x) unlist(colnames(interactions[,12:94])[interactions[x,12:94] == 1]))
names(cell.interactions) <- interactions$interacting_pair

names(cell.interactions)


sigmeans <- read.delim('cpdb/degs_analysis/significant_means.txt')
unlist(strsplit(sigmeans$interacting_pair, '_')) %in% rownames(chrX)
sigmeans[7,]