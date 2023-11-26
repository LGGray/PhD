library(edgeR)
library(Seurat)
library(qvalue)
library(tidyverse)
library(ggplot2)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(reshape2)
library(clusterProfiler)

if(dir.exists('differential.expression/sex_interaction') != TRUE){dir.create('differential.expression/sex_interaction')}

# Read in file from command line
pbmc <- readRDS(commandArgs(trailingOnly = TRUE)[1])
assay <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

# if('sex' %in% colnames(pbmc@meta.data) == FALSE){
#     print('sex not in metadata...predicting sex')
#     source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/predict.sex.R')
#     pbmc <- predict.sex(pbmc, assay='decontXcounts', slot='data', individual='individual')
# }

pbmc <- subset(pbmc, disease_state %in% c('na', 'managed'))
pbmc$age <- as.numeric(gsub('-year-old human stage', '', pbmc$development_stage))

for (cell in levels(pbmc)){
  # Select cell type
  print(cell)
  # subset object by cell type
  pbmc.cell <- subset(pbmc, cellTypist == cell)
  # check if there are enough cells and skip if not
  if(nrow(pbmc.cell) < 10){
    print("Not enough cells")
    next
  }
  # check if there are enough cell in both conditions and skip if not
  if(length(unique(pbmc.cell$condition)) != 2){
    print("Not enough conditions")
    next
  } else {
    # print(table(pbmc.cell$condition, pbmc.cell$cellTypist))
  }
  # Keep genes with expression in 5% of cells
  keep <- rowSums(pbmc.cell@assays$RNA@counts > 0) > ncol(pbmc.cell) * 0.05
  features <- names(keep[keep == T])
  pbmc.cell <- subset(pbmc.cell, features=features)

  # Calculate cellcount
  cellCount <- as.data.frame.matrix(table(pbmc.cell$individual, pbmc.cell$cellTypist))

  # Psudobulking by summing counts
  expr <- AggregateExpression(pbmc.cell, group.by='individual', slot='counts')[[assay]]
  expr <- expr[(rowSums(expr) > 0),]

  # edgeR-QLFTest
  targets = unique(data.frame(condition = pbmc.cell$condition,
                      individual = pbmc.cell$individual,
                      sex = pbmc.cell$sex,
                      age = pbmc.cell$age))
  targets$cellCount <- cellCount[,cell]
  targets <- targets[match(colnames(expr), targets$individual),]
  rownames(targets) <- targets$individual
  targets$sex <- factor(targets$sex)
  targets$sex <- relevel(targets$sex, ref = "F")
  targets$condition <- factor(targets$condition)
  targets$condition <- relevel(targets$condition, ref = "disease")
  design <- model.matrix(~0 + cellCount + sex + condition + age + sex:condition, data=targets)
  y = DGEList(counts = expr, group = targets$condition)
  # y <- calcNormFactors(y)
  y = estimateGLMRobustDisp(y, design, trend.method = 'auto')
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef=6)
  print(summary(decideTests(qlf)))
  res = topTags(qlf, n = Inf) %>%
    as.data.frame() %>%
    rownames_to_column('gene')
  res$FDR <- qvalue(p = res$PValue)$qvalues
  cell = gsub("/|-| ", "_", cell)
  write.table(res, paste0("differential.expression/age_sex_interaction/", cell, ".txt"),
              row.names=F, sep="\t", quote = F)
}

meta <- unique(pbmc@meta.data[, c('individual', 'condition', 'sex', 'age')])
table(meta$condition, meta$sex)


load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/fishers.test.degs.R')

condition <- deg.list('differential.expression/edgeR_age', logfc=0.2)
names(condition) <- replace.names(gsub('_', '.', names(condition)))
sex <- deg.list('differential.expression/age_sex_interaction', logfc=0.2)
names(sex) <- replace.names(gsub('_', '.', names(sex)))

intersection <- lapply(names(sex), function(x){
  merge(sex[[x]], condition[[x]], by='gene', all=T)
})
names(intersection) <- names(sex)

intersection.chrX <- lapply(intersection, function(x){
  subset(x, gene %in% rownames(chrX))
})

lapply(intersection.chrX, function(x){
  tmp <- subset(x, !is.na(logFC) & !is.na(logFC.disease_vs_control))
  cor(tmp$logFC, tmp$logFC.disease_vs_control, )
})

lapply(sex, nrow)


fisher.test.edgeR(res, rownames(chrX), logfc=0.5, direction='none')

res[abs(res$logFC) > 0.5 & res$FDR < 0.05,'gene'] %in% rownames(chrX)

