library(edgeR)
library(Seurat)
library(qvalue)

setwd("~/datasets/integrated/")

args = commandArgs(trailingOnly=TRUE)
args = as.numeric(args[1])

pbmc <- readRDS("immune.combined.RDS")
# Add overall disease and control
pbmc$condition <- ifelse(pbmc$disease %in% c('Control', 'OA'), 'Control', 'Disease')
# Select cell type
cell = levels(pbmc)[args]
print(cell)
# subset object by cell type
pbmc.cell <- subset(pbmc, predicted.celltype.l2 == cell)
expr <- GetAssayData(object = pbmc.cell, assay='SCT', slot = "counts")
expr = as.data.frame(expr)

# edgeR-LRT
targets = data.frame(group = pbmc.cell$condition,
                     individual = pbmc.cell$individual)
design <- model.matrix(~group, data=targets)
y = DGEList(counts = expr, group = targets$group)
# Disease group as reference
y$samples$group <- factor(y$samples$group, levels = c('Disease', 'Control'))
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
write.table(lrt, paste0("edgeR-LRT/", cell, ".txt"),
            row.names=F, sep="\t", quote = F)
