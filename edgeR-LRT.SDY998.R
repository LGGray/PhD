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
expr <- GetAssayData(object = pbmc.cell, slot = "counts")
expr = as.data.frame(expr)

# edgeR-LRT
targets = data.frame(group = pbmc.cell$disease,
                     individual = pbmc.cell$sample,
                     lane = pbmc.cell$lane)
design <- model.matrix(~lane+group, data=targets)
y = DGEList(counts = expr, group = targets$group)
# Disease group as reference
y$samples$group <- factor(y$samples$group, levels = c('RA', 'OA'))
# Filter for expression in 5% of cells
keep <- rowSums(expr > 0) > dim(y)[2]*0.05
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design)
lrt <- glmLRT(fit)
print(summary(decideTests(lrt)))
lrt <- as.data.frame(lrt)
lrt$FDR <- qvalue(p = lrt$PValue)$qvalues
gene <- rownames(lrt)
lrt <- cbind(gene,lrt)
cell = sub(" ", "_", cell)
write.table(lrt, paste0("~/datasets/SDY998/edgeR-LRT/", cell, ".txt"),
            row.names=F, sep="\t", quote = F)