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
batch <- commandArgs(trailingOnly = TRUE)[2]

for (cell in levels(pbmc)){
  # Select cell type
  print(cell)
  # subset object by cell type
  pbmc.cell <- subset(pbmc, cellTypist == cell)
  # check if there are enough cells and skip if not
  if(nrow(pbmc.cell) < 30){
    print("Not enough cells")
    next
  }
  # check if there are enough cell in both conditions and skip if not
  if(length(unique(pbmc.cell$condition)) != 2){
    print("Not enough conditions")
    next
  } else {
    table(pbmc.cell$condition, pbmc.cell$cellTypist)
  }
  # Keep genes with expression in 5% of cells
  keep <- rowSums(pbmc.cell@assays$RNA@counts > 0) > ncol(pbmc.cell) * 0.05
  features <- names(keep[keep == T])
  pbmc.cell <- subset(pbmc.cell, features=features)

  # Psudobulking by summing counts
  expr <- AggregateExpression(pbmc.cell, group.by='individual', slot='counts')$RNA
  expr <- expr[,(colSums(expr) > 0)]
  # edgeR-QLFTest
  targets = unique(data.frame(group = pbmc.cell$condition,
                      individual = pbmc.cell$individual,
                      batch = pbmc.cell@meta.data[,batch]))
  targets <- targets[match(colnames(expr), targets$individual),]
  # targets$SV1 <- tapply(pbmc.cell$SV1, pbmc.cell$individual, sum)
  # targets$SV2 <- tapply(pbmc.cell$SV2, pbmc.cell$individual, sum)
  rownames(targets) <- targets$individual
  design <- model.matrix(~0 + batch + group, data=targets)
  y = DGEList(counts = expr, group = targets$group)
  # Disease group as reference
  contrasts <- makeContrasts(disease_vs_control = groupdisease - groupcontrol, levels = design)
  #y <- calcNormFactors(y, method='TMM')
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

#  MAST
for (cell in levels(pbmc)){
  # Select cell type
  print(cell)
  # subset object by cell type
  pbmc.cell <- subset(pbmc, cellTypist == cell)
  # Keep genes with expression in 5% of cells
  keep <- rowSums(pbmc.cell@assays$RNA@counts > 0) > ncol(pbmc.cell) * 0.05
  features <- names(keep[keep == T])
  pbmc.cell <- subset(pbmc.cell, features=features)

  # check if there are enough cell in both conditions and skip if not
  if(length(unique(pbmc.cell$condition)) != 2){
    print("Not enough conditions")
    next
  } else {
    table(pbmc.cell$condition, pbmc.cell$cellTypist)
  }
  Idents(pbmc.cell) <- "condition"
  result <- FindMarkers(pbmc.cell, slot='counts', ident.1 = "disease", ident.2 = "control",
                             test.use = "MAST", latent.vars=batch, min.pct = 0, logfc.threshold = 0)                          
  result <- cbind(gene = rownames(result), result)
  cell = gsub("/|-| ", "_", cell)
  write.table(result, paste0("differential.expression/MAST/", cell, ".txt"),
              row.names=F, sep="\t", quote = F)
}

print("Done with MAST")