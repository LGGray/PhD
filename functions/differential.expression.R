library(edgeR)
library(Seurat)
library(MAST)
library(qvalue)
library(tidyverse)

if(dir.exists('psuedobulk') != TRUE){dir.create('psuedobulk')}

pbmc <- readRDS("pbmc.female.RDS")
metadata <- 
batch <- split(pbmc@meta.data, pbmc@meta.data$individual) %>%
  lapply(function(x) data.frame(batch_1=sum(x$batch_1), batch_2=sum(x$batch_2))) %>%
  bind_rows(.id='individual')


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

  # Psuedobulk
  individual <- as.factor(pbmc.cell$individual)
  mm <- model.matrix(~ 0 + individual)
  colnames(mm) <- levels(individual)
  expr <- GetAssayData(object = pbmc.cell, slot = "counts") %*% mm

  # edgeR-QLFTest
  targets = unique(data.frame(group = pbmc.cell$condition,
                      individual = pbmc.cell$individual))
  targets <- targets[match(colnames(expr), targets$individual),]
  rownames(targets) <- targets$individual
  targets <- merge(targets, batch, by='individual')
  design <- model.matrix(~0+ group, data=targets)
  y = DGEList(counts = expr, group = targets$group)
  # Disease group as reference
  contrasts <- makeContrasts(disease_vs_control = groupdisease - groupcontrol, levels = design)
  #y <- calcNormFactors(y, method='TMM')
  y = estimateGLMRobustDisp(y, design,
                            trend.method = 'auto')
  fit <- glmQLFit(y, design)
  contrast_matrix <- contrasts[ ,c("disease_vs_control")]
  qlf <- glmQLFTest(fit, contrast=contrast_matrix)
  print(summary(decideTests(qlf)))
  res = topTags(qlf, n = Inf) %>%
    as.data.frame() %>%
    rownames_to_column('gene')
  res$FDR <- qvalue(p = res$PValue)$qvalues
  cell = gsub("/|-| ", "_", cell)
  write.table(res, paste0("psuedobulk/", cell, ".edgeR-QLF.txt"),
              row.names=F, sep="\t", quote = F)
}

print("Done with edgeR-QLF")

# #  Seurat Wilcoxon rank sum test
# for (cell in levels(pbmc)){
#   # Select cell type
#   print(cell)
#   # subset object by cell type
#   pbmc.cell <- subset(pbmc, cellTypist == cell)

#   # check if there are enough cell in both conditions and skip if not
#   if(length(unique(pbmc.cell$condition)) != 2){
#     print("Not enough conditions")
#     next
#   } else {
#     table(pbmc.cell$condition, pbmc.cell$cellTypist)
#   }
#   Idents(pbmc.cell) <- "condition"
#   result <- FindMarkers(pbmc.cell, slot='counts', ident.1 = "disease", ident.2 = "control",
#                              test.use = "wilcox", min.pct = 0, logfc.threshold = 0)
#   result <- cbind(gene = rownames(result), result)

#   write.table(result, paste0("psuedobulk/", cell, ".wilcoxon.txt"),
#               row.names=F, sep="\t", quote = F))
# }

#  MAST
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
  Idents(pbmc.cell) <- "condition"
  result <- FindMarkers(pbmc.cell, slot='counts', ident.1 = "disease", ident.2 = "control",
                             test.use = "MAST", latent.vars=c('batch_1', 'batch_2'), min.pct = 0, logfc.threshold = 0)                          
  result <- cbind(gene = rownames(result), result)
  cell = gsub("/|-| ", "_", cell)
  write.table(result, paste0("psuedobulk/", cell, ".MAST.txt"),
              row.names=F, sep="\t", quote = F)

print("Done with MAST")

# mast <- result
# qlf_wilcox <- merge(qlf, wilcox, by = "gene", all = T)
# qlf_mast <- merge(qlf, mast, by = "gene", all = T)
# wilcox_mast <- merge(wilcox, mast, by = "gene", all = T)

# cor(qlf_wilcox$logFC, qlf_wilcox$avg_log2FC)
# cor(qlf_mast$logFC, qlf_mast$avg_log2FC)
# cor(wilcox_mast$avg_log2FC.x, wilcox_mast$avg_log2FC.y)