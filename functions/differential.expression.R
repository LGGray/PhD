library(edgeR)
library(Seurat)
library(MAST)
library(qvalue)
library(tidyverse)

if(dir.exists('differential.expression') != TRUE){dir.create('differential.expression')}
if(dir.exists('differential.expression/edgeR') != TRUE){dir.create('differential.expression/edgeR')}
if(dir.exists('differential.expression/MAST') != TRUE){dir.create('differential.expression/MAST')}

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
  targets <- targets[match(colnames(expr), targets$individual),]
  rownames(targets) <- targets$individual
  design <- model.matrix(~0 + group, data=targets)
  y = DGEList(counts = expr, group = targets$group)
  # Disease group as reference
  contrasts <- makeContrasts(disease_vs_control = groupdisease - groupcontrol, levels = design)
  y <- calcNormFactors(y, method='TMM')
  y = estimateGLMRobustDisp(y, design, trend.method = 'auto')
  fit <- glmQLFit(y, design)
  contrast_matrix <- contrasts[ ,c("disease_vs_control")]
  qlf <- glmQLFTest(fit, contrast=contrast_matrix)
  print(summary(decideTests(qlf)))
  res = topTags(qlf, n = Inf) %>%
    as.data.frame() %>%
    rownames_to_column('gene')
  res$FDR <- qvalue(p = res$PValue)$qvalues
  cell = gsub("/|-| ", "_", cell)
  write.table(res, paste0("differential.expression/edgeR/", cell, ".txt"),
              row.names=F, sep="\t", quote = F)
}

print("Done with edgeR-QLF")