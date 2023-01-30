library(Seurat)
library(ggplot2)
library(dplyr)

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/chisq.test.degs.R')

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/RA_SDY998/')

# Load the data
pbmc <- readRDS('pbmc.female.RDS')

deg.list <- edgeR.list('psuedobulk', logfc=0.5)
names(deg.list) <- gsub('.edgeR-LRT', '', names(deg.list))

deg.metrics <- lapply(deg.list, function(x){
    up.all <- x[which(x$logFC > 0),]
    down.all <- x[which(x$logFC < 0),]
    XCI <- subset(x, gene %in% rownames(chrX))
    up.XCI <- XCI[which(XCI$logFC > 0),]
    down.XCI <- XCI[which(XCI$logFC < 0),]
    data.frame(up=nrow(up.all), down=nrow(down.all), up.chrX=nrow(up.XCI), down.chrX=nrow(down.XCI))
}) %>% do.call(rbind, .) 
deg.metrics

# Celltype abundance
levels(pbmc)
table(pbmc$cellTypist, pbmc$condition)


