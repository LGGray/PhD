library(Seurat)
library(dplyr)
library(tibble)
library(randomForest)
library(Boruta)
library(caret)
library(pROC)
library(ggplot2)
source('~/R_code/functions/edgeR.list.R')
load('~/datasets/XCI/escapees.Rdata')

setwd('~/datasets/integrated/')
pbmc <- readRDS('immune.combined.RDS')
edgeR <- edgeR.list('psuedobulk/', logfc=0.5)

args = commandArgs(trailingOnly=TRUE)
args <- as.numeric(args)
cell <- levels(pbmc)[args]

xcape <- subset(edgeR[[gsub(' ', '_', cell)]], gene %in% rownames(escape))
features <- xcape$gene

pbmc.subset <- subset(pbmc, predicted.celltype.l2 %in% cell)

expr <- t(as.matrix(GetAssayData(pbmc.subset, assay='SCT')))
expr <- data.frame(expr[,which(colnames(expr) %in% features)])
expr$condition <- factor(pbmc.subset$condition, levels = c('Control', 'Disease'))

set.seed(42)
# Split data into training and test sets 80:20
training <- expr %>%
  rownames_to_column('cell') %>%
  group_by(condition) %>%
  sample_frac(0.8) %>%
  column_to_rownames('cell')
training.idx <- rownames(expr) %in% rownames(training)

testing <- expr[-training.idx,]

# Feature Selection
boruta <- Boruta(condition ~ ., data = expr, doTrace = 2, maxRuns = 500)
final.boruta <- TentativeRoughFix(boruta)
# Build randomForest model
rf <- randomForest(
  condition ~ .,
  data=training,
  importance=T
)

pred <- predict(rf, newdata=testing)
cm <- confusionMatrix(testing$condition, pred)
print(cm)

new.formula <- getNonRejectedFormula(final.boruta)
rfboruta <- randomForest(
  new.formula,
  data=training,
  importance=T
)
save(rfboruta, file= paste0('feature.selection/randomForest/', sub(' ', '_', cell), '.model.RData'))
new.pred <- predict(rfboruta, newdata=testing)
new.cm <- confusionMatrix(testing$condition, new.pred)
print(new.cm)

new.features <- names(final.boruta$finalDecision[which(final.boruta$finalDecision == 'Confirmed')])
df <- subset(xcape, gene %in% new.features)
write.table(df, paste0('feature.selection/randomForest/', sub(' ', '_', cell), '.txt'), row.names=F, quote=F, sep='\t')

# plot ROC curve
predicted <- apply(predict(rfboruta, newdata=testing, type='prob'), 1, max)
rocobj <- roc(testing$condition, predicted)
auc <- round(auc(testing$condition, predicted),4)
pdf(paste0('feature.selection/randomForest/', sub(' ', '_', cell), '.ROC.pdf'))
ggroc(rocobj, size=2) +
  ggtitle(paste0(cell, ': ROC Curve ', '(AUC = ', auc, ')')) +
  theme_minimal()
dev.off()
