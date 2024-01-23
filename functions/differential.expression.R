library(edgeR)
library(Seurat)
library(MAST)
library(qvalue)
library(tidyverse)
library(ggplot2)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(reshape2)
library(clusterProfiler)

if(dir.exists('differential.expression') != TRUE){dir.create('differential.expression')}
if(dir.exists('differential.expression/edgeR_cellCount') != TRUE){dir.create('differential.expression/edgeR_cellCount')}

# Read in file from command line
pbmc <- readRDS(commandArgs(trailingOnly = TRUE)[1])

for (cell in levels(pbmc)){
  # Select cell type
  print(cell)
  # subset object by cell type
  pbmc.cell <- subset(pbmc, cellTypist == cell)

  # Remove genes with pseudobulked expression in less than 5% of individuals
  expr <- AverageExpression(pbmc.cell, group.by='individual', slot='counts')$RNA
  keep <- apply(expr, 1, function(x) sum(x > 0) > ncol(expr) * 0.05)
  expr <- expr[keep,]

  targets = unique(data.frame(group = pbmc.cell$condition,
                      individual = pbmc.cell$individual))
  targets$cellCount <- sapply(targets$individual, function(id) sum(pbmc.cell$individual == id))
  targets <- targets[match(colnames(expr), targets$individual),]
  rownames(targets) <- targets$individual
  design <- model.matrix(~0 + cellCount + group, data=targets)
  y = DGEList(counts = expr, group = targets$group)
  # Disease group as reference
  contrasts <- makeContrasts(disease_vs_control = groupdisease - groupcontrol,
                            cell_type_effect = cellCount, levels = design)
  y <- calcNormFactors(y, method='TMM')
  y = estimateGLMRobustDisp(y, design, trend.method = 'auto')
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, contrast=contrasts)
  print(summary(decideTests(qlf)))
  res = topTags(qlf, n = Inf) %>%
    as.data.frame() %>%
    rownames_to_column('gene')
  res$FDR <- qvalue(p = res$PValue)$qvalues
  cell = gsub("/|-| ", "_", cell)
  write.table(res, paste0("differential.expression/edgeR_cellCount/", cell, ".txt"),
              row.names=F, sep="\t", quote = F)
}

print("Done with edgeR-QLF")