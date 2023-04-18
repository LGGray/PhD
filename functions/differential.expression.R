library(edgeR)
library(Seurat)
library(qvalue)
library(tidyverse)

if(dir.exists('psuedobulk') != TRUE){dir.create('psuedobulk')}

pbmc <- readRDS("pbmc.female.RDS")

# Keep genes with expression in 5% of cells
exp <- GetAssayData(pbmc, slot = "counts")
keep <- apply(exp, 1, function(x) sum(x > 0) > ncol(pbmc)*0.05)
features <- names(keep[keep == T])
pbmc <- subset(pbmc, features=features)

for (cell in levels(pbmc)){
  # Select cell type
  print(cell)
  # subset object by cell type
  pbmc.cell <- subset(pbmc, cellTypist == cell)

  # check if there are enough cell in both conditions and skip if not
  if(length(unique(pbmc.cell$condition)) != 2){
    print("Not enough conditions")
    next
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
  lrt <- glmQLFTest(fit)
  print(summary(decideTests(lrt)))
  res = topTags(lrt, n = Inf) %>%
    as.data.frame() %>%
    rownames_to_column('gene')
  res$FDR <- qvalue(p = res$PValue)$qvalues
  cell = gsub("/|-| ", "_", cell)
  write.table(res, paste0("psuedobulk/", cell, ".edgeR-QLF.txt"),
              row.names=F, sep="\t", quote = F)
}