library(Seurat)
library(edgeR)
library(qvalue)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
args = as.numeric(args[1])

pbmc <- readRDS('/directflow/SCCGGroupShare/projects/lacgra/lupus.Chun/asian.female.RDS')
# Select cell type
cell = levels(pbmc)[args]
print(cell)
# subset object by cell type
pbmc.cell <- subset(pbmc, predicted.celltype.l2 == cell)

pbmc.cell$age <- as.numeric(gsub('-year-old human stage', '', pbmc.cell$development_stage))

# Psuedobulk
sample <- as.factor(pbmc.cell$ind_cov)
mm <- model.matrix(~0+sample)
colnames(mm) <- levels(sample)
expr <- GetAssayData(object = pbmc.cell, slot = "counts") %*% mm

# edgeR-LRT
targets = unique(data.frame(group = pbmc.cell$disease,
                            age = pbmc.cell$age,
                            pool = pbmc.cell$Processing_Cohort,
                            individual = pbmc.cell$ind_cov))
targets <- targets[match(colnames(expr), targets$individual),]
design <- model.matrix(~age+pool+group, data=targets)
y = DGEList(counts = expr, group = targets$group)
# Disease group as reference
y$samples$group <- relevel(y$samples$group, ref='normal')
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
write.table(res, paste0("/directflow/SCCGGroupShare/projects/lacgra/lupus.Chun/pseudobulk/asian/", cell, ".txt"),
            row.names=F, sep="\t", quote = F)
