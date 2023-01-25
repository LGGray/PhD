library(edgeR)
library(Seurat)
library(qvalue)
library(tidyverse)

setwd("/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/RA_SDY998")

if(dir.exists('psuedobulk') != TRUE){dir.create('psuedobulk')}

pbmc <- readRDS("pbmc.female.RDS")
DefaultAssay(pbmc) <- 'SCT' 

for (cell in levels(pbmc)){
  # Select cell type
  print(cell)
  # subset object by cell type
  pbmc.cell <- subset(pbmc, cellTypist == cell)

  # check if there are enough cell in both conditions and stop if not
  if(length(unique(pbmc.cell$condition)) != 2){
    print("Not enough conditions")
    stop()
  } else {
    table(pbmc.cell$condition, pbmc.cell$cellTypist)
  }

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
  y$samples$group <- factor(y$samples$group, levels = c('disease', 'control'))
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
  write.table(res, paste0("psuedobulk/", cell, ".edgeR-LRT.txt"),
              row.names=F, sep="\t", quote = F)
}