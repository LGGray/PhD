library(edgeR)
library(Seurat)
library(qvalue)

setwd("~/datasets/SDY998/")

args = commandArgs(trailingOnly=TRUE)
args = as.numeric(args)

pbmc <- readRDS("pbmc.female.RDS")
# Select cell type
cell = levels(pbmc)[args]
print(cell)
# subset object by cell type
pbmc.cell <- subset(pbmc, predicted.celltype.l2 == cell)

# Psuedobulk
sample <- as.factor(pbmc.cell$sample)
mm <- model.matrix(~ 0 + sample)
colnames(mm) <- levels(sample)
expr <- GetAssayData(object = pbmc.cell, slot = "counts") %*% mm


# edgeR-LRT
targets = unique(data.frame(group = pbmc.cell$disease,
                     individual = pbmc.cell$sample))
design <- model.matrix(~0+group, data=targets)
y = DGEList(counts = expr, group = targets$group)
# Disease group as reference
y$samples$group <- factor(y$samples$group, levels = c('RA', 'OA'))
# Filter for expression in 5% of cells
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design)
lrt <- glmLRT(fit)
print(summary(decideTests(lrt)))
lrt <- as.data.frame(lrt$table)
lrt$FDR <- p.adjust(lrt$PValue, method = 'fdr')
gene <- rownames(lrt)
lrt <- cbind(gene,lrt)
cell = sub(" ", "_", cell)
write.table(lrt, paste0("~/datasets/SDY998/psuedobulk/", cell, ".txt"),
            row.names=F, sep="\t", quote = F)
