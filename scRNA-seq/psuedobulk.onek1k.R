library(Seurat)
library(edgeR)
library(qvalue)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
args = as.numeric(args[1])

pbmc <- readRDS('/directflow/SCCGGroupShare/projects/lacgra/seurat.object/onek1k.RDS')
# Select cell type
cell = levels(pbmc)[args]
print(cell)
# subset object by cell type
pbmc.cell <- subset(pbmc, predicted.celltype.l3 == cell)

# Psuedobulk
sample <- as.factor(pbmc.cell$individual)
mm <- model.matrix(~0+sample)
colnames(mm) <- levels(sample)
expr <- GetAssayData(object = pbmc.cell, slot = "counts") %*% mm

# edgeR-LRT
targets = unique(data.frame(group = pbmc.cell$sex,
                            age = pbmc.cell$age,
                            pool = pbmc.cell$pool,
                            individual = pbmc.cell$individual))
targets <- targets[match(colnames(expr), targets$individual),]
design <- model.matrix(~age+group, data=targets)
y = DGEList(counts = expr, group = targets$group)
# Disease group as reference
y$samples$group <- relevel(y$samples$group, ref='F')
# Filter for expression in 5% of cells
y <- calcNormFactors(y, method='TMM')
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design)
lrt <- glmLRT(fit)
print(summary(decideTests(lrt)))
res = topTags(lrt, n = Inf) %>%
  as.data.frame() %>%
  rownames_to_column('gene')
res$FDR <- qvalue(p = res$PValue)$qvalues
cell = sub(" ", "_", cell)
write.table(res, paste0("/directflow/SCCGGroupShare/projects/lacgra/sex.bias/edgeR-LRT/", cell, ".txt"),
            row.names=F, sep="\t", quote = F)
