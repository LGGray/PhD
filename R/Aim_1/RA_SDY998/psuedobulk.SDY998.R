library(edgeR)
library(Seurat)
library(qvalue)
library(tidyverse)

setwd("/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/RA_SDY998")

if(dir.exists('psuedobulk') != TRUE){dir.create('psuedobulk')}

args = commandArgs(trailingOnly=TRUE)
args = as.numeric(args)

pbmc <- readRDS("pbmc.female.RDS")
DefaultAssay(pbmc) <- 'SCT' 

# Select cell type
cell = levels(pbmc)[args]
print(cell)
# subset object by cell type
pbmc.cell <- subset(pbmc, predicted.celltype.l2 == cell)

# Psuedobulk
individual <- as.factor(pbmc.cell$individual)
mm <- model.matrix(~ 0 + individual)
colnames(mm) <- levels(individual)
expr <- GetAssayData(object = pbmc.cell, slot = "counts") %*% mm

# edgeR-LRT
targets = unique(data.frame(group = pbmc.cell$condition,
                     individual = pbmc.cell$individual))
targets <- targets[match(colnames(expr), targets$individual),]
design <- model.matrix(~group, data=targets)
y = DGEList(counts = expr, group = targets$group)
# Disease group as reference
y$samples$group <- factor(y$samples$group, levels = c('RA', 'OA'))
# Filter for expression in 5% of cells
y <- calcNormFactors(y, method='TMM')
y = estimateGLMRobustDisp(y, design,
                          trend.method = 'auto')
fit <- glmQLFit(y, design)
lrt <- glmLRT(fit)
print(summary(decideTests(lrt)))
res = topTags(lrt, n = Inf) %>%
  as.data.frame() %>%
  rownames_to_column('gene')
res$FDR <- qvalue(p = res$PValue)$qvalues
cell = sub(" ", "_", cell)
write.table(res, paste0("psuedobulk/", cell, ".txt"),
            row.names=F, sep="\t", quote = F)
