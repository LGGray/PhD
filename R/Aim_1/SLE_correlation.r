library(Seurat)
library(speckle)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(cluster)
library(ggrepel)
library(fgsea)
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/fishers.test.degs.R')

# Load colours 
load('/directflow/SCCGGroupShare/projects/lacgra/PhD/R/celltype.colours.RData')
load('/directflow/SCCGGroupShare/projects/lacgra/PhD/R/study_colours.RData')

# Read in gene sets 
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')
X.immune <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/X.immune.txt')[,1]
chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- chrX$Gene.name

# Read in edgeR files for each study
MS_all <- deg.list('MS_GSE193770/differential.expression/edgeR', filter=FALSE)
pSS_all <- deg.list('pSS_GSE157278/differential.expression/edgeR', filter=FALSE)
UC_all <- deg.list('UC_GSE125527/differential.expression/edgeR', filter=FALSE)
CO_all <- deg.list('CD_Kong/colon/differential.expression/edgeR', filter=FALSE)
names(CO_all)[2] <- 'CD16-_NK_cells'
TI_all <- deg.list('CD_Kong/TI/differential.expression/edgeR', filter=FALSE)
names(TI_all)[1] <- 'CD16-_NK_cells'
SLE_all <- deg.list('lupus_Chun/differential.expression/edgeR', filter=FALSE)


common_celltypes <- Reduce(intersect, list(names(pSS_all), names(UC_all), names(CO_all), names(TI_all), names(SLE_all)))

results_list <- list()
for(cell in common_celltypes){
    SLE_deg <- subset(SLE_all[[cell]], FDR < 0.05 & abs(logFC) > 0.1)

    pSS <- merge(SLE_deg, pSS_all[[cell]], by='gene', suffixes=c('_SLE', '_pSS'))
    pSS_correlation <- cor.test(pSS$logFC_SLE, pSS$logFC_pSS, method='spearman')
    UC <- merge(SLE_deg, UC_all[[cell]], by='gene', suffixes=c('_SLE', '_UC'))
    UC_correlation <- cor.test(UC$logFC_SLE, UC$logFC_UC, method='spearman')
    CO <- merge(SLE_deg, CO_all[[cell]], by='gene', suffixes=c('_SLE', '_CO'))
    CO_correlation <- cor.test(CO$logFC_SLE, CO$logFC_CO, method='spearman')
    TI <- merge(SLE_deg, TI_all[[cell]], by='gene', suffixes=c('_SLE', '_TI'))
    TI_correlation <- cor.test(TI$logFC_SLE, TI$logFC_TI, method='spearman')

    results <- data.frame(cell=cell, pSS_correlation=pSS_correlation$estimate, 
    pSS_p=pSS_correlation$p.value, UC_correlation=UC_correlation$estimate, 
    UC_p=UC_correlation$p.value, CO_correlation=CO_correlation$estimate, 
    CO_p=CO_correlation$p.value, TI_correlation=TI_correlation$estimate, 
    TI_p=TI_correlation$p.value, row.names='')
    results_list[[cell]] <- results
}

results <- do.call(rbind, results_list)
write.csv(results, 'SLE_study_correlation.csv', row.names=FALSE)


