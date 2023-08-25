library(edgeR)
library(Seurat)
library(qvalue)
library(tibble)

if(dir.exists('differential.expression') != TRUE){dir.create('differential.expression')}
if(dir.exists('differential.expression/edgeR') != TRUE){dir.create('differential.expression/edgeR')}
if(dir.exists('differential.expression/MAST') != TRUE){dir.create('differential.expression/MAST')}

# Read in Seurat file
pbmc <- readRDS('pbmc.female.control-managed.RDS')
pbmc$disease_state <- gsub('^na', 'control', pbmc$disease_state)
pbmc$development_stage <- as.numeric(gsub('-.+', '', pbmc$development_stage))

cell <- levels(pbmc)[as.numeric(commandArgs(trailingOnly = TRUE)[1])]

# Select cell type
print(cell)
# subset object by cell type
pbmc.cell <- subset(pbmc, cellTypist == cell)
rm(pbmc)
# Keep genes with expression in 5% of cells
keep <- rowSums(pbmc.cell@assays$RNA@counts > 0) > ncol(pbmc.cell) * 0.05
features <- names(keep[keep == T])
pbmc.cell <- subset(pbmc.cell, features=features)

# Psudobulking by summing counts
expr <- AggregateExpression(pbmc.cell, group.by='individual', slot='counts')$RNA
expr <- expr[,(colSums(expr) > 0)]
# edgeR-QLFTest
abundance <- as.data.frame.matrix(table(pbmc.subset$individual, pbmc.subset$cellTypist), row.names=NULL)
colnames(abundance) <- 'celltype'
groups <- cbind(groups, abundance)
targets = unique(data.frame(individual = pbmc.cell$individual,
                            group = pbmc.cell$condition,
                            age = pbmc.cell$development_stage))
targets <- cbind(targets, abundance)
rownames(targets) <- NULL
targets <- targets[match(colnames(expr), targets$individual),]
design <- model.matrix(~ 0 + age + celltype + condition, data=targets)
contrasts <- makeContrasts(disease_vs_control = conditiondisease - conditioncontrol,
                           age_effect = development_stage,
                           cell_type_effect = celltype,
                           levels = design) 
y = DGEList(counts = expr, group = targets$group)
# Disease group as reference
contrasts <- makeContrasts(disease_vs_control = groupdisease - groupcontrol, levels = design)
y = estimateGLMRobustDisp(y, design, trend.method = 'auto')
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, contrast=contrasts)
print(summary(decideTests(qlf)))
res = topTags(qlf, n = Inf) %>%
    as.data.frame() %>%
    rownames_to_column('gene')
    res$FDR <- qvalue(p = res$PValue)$qvalues
cell = gsub("/|-| ", "_", cell)
write.table(res, paste0("differential.expression/edgeR/", cell, ".txt"),
            row.names=F, sep="\t", quote = F)