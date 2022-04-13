library(edgeR)
library(Seurat)
library(qvalue)

args = commandArgs(trailingOnly=TRUE)
args = as.numeric(args[1])

pbmc <- readRDS("/directflow/SCCGGroupShare/projects/lacgra/seurat.object/pbmc.all.RDS")
# Select cell type
cell = levels(pbmc)[args]
print(cell)
# subset object by cell type
pbmc.cell <- subset(pbmc, predicted.celltype.l2 %in% cell)
rm(pbmc)
expr <- GetAssayData(object = pbmc.cell, slot = "counts")
expr = as.data.frame(expr)
pbmc.cell$individual <- factor(pbmc.cell$individual)
# edgeR-LRT
targets = data.frame(group = as.factor(pbmc.cell$RA, levels(c('Y', 'N')),
                     pool = pbmc.cell$pool,
                     age = pbmc.cell$age,
                     individual = pbmc.cell$individual)
design <- model.matrix(~0+pool+age+group, data=targets)
y = DGEList(counts = expr, group = targets$group)
                     
# Filter for expression in 5% of cells
keep <- rowSums(expr > 0) > dim(y)[2]*0.05
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design)
lrt <- glmLRT(fit, coef=dim(design)[2])
print(summary(decideTests(lrt)))
lrt <- as.data.frame(lrt)
lrt$FDR <- qvalue(p = lrt$PValue)$qvalues
gene <- rownames(lrt)
lrt <- cbind(gene,lrt)
cell = sub(" ", "_", cell)
write.table(lrt, paste0("~/datasets/OneK1k/C_vs_RA_DEout/edgeR-LRT/all.sc/", cell, ".txt"),
            row.names=F, sep="\t", quote = F)
