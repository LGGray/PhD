library(BDA)
library(SeuratObject)
library(Seurat)

pbmc <- readRDS(commandArgs(trailingOnly = TRUE)[1])

results <- BDA(object = pbmc,
               contrast.by = "condition", ident.1 = "disease", ident.2 = "control",
               group.by = "cellTypist", cell.populations = levels(pbmc), 
               p.adjust.method = "fdr", min.pct = 0.05, n.cores = 4, method = "LR")

save(results, file = "BDA.RDS")