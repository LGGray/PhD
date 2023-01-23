library(Seurat)
library(dplyr)
library(tibble)
library(MASS)
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
training.data <- expr %>%
  rownames_to_column('cell') %>%
  group_by(condition) %>%
  sample_frac(0.8) %>%
  column_to_rownames('cell')
training <- rownames(expr) %in% rownames(training.data)

testing <- expr[-training,]

# Fit full model
full.model <- glm.fits <- glm(
  condition ~ .,
  data = expr, family = binomial, subset = training
)
coef(full.model)
step.model <- full.model %>% stepAIC(trace = F)
coef(step.model)
# save model to file
save(step.model, file= paste0('feature.selection/logit/', sub(' ', '_', cell), '.model.RData'))

# Prediction accuracy of full model
probabilities <- full.model %>% predict(testing, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "Disease", "Control")
print('Full model accuracy: ')
print(mean(predicted.classes == testing$condition))

# Prediction accuracy of stepwise model
probabilities <- predict(step.model, testing, type='response')
predicted.classes <- ifelse(probabilities > 0.5, "Disease", "Control")
print('stepwise model accuracy: ')
print(mean(predicted.classes == testing$condition))

# save genes to file
df <- subset(xcape, gene %in% names(coef(step.model)))
write.table(df, paste0('feature.selection/', sub(' ', '_', cell), '.txt'), row.names=F, quote=F, sep='\t')

# plot ROC curve
predicted <- predict(step.model, newdata=testing, type='response')
rocobj <- roc(testing$condition, predicted)
auc <- round(auc(testing$condition, predicted),4)
pdf(paste0('feature.selection/logit/', sub(' ', '_', cell), '.ROC.pdf'))
ggroc(rocobj, size=2) +
  ggtitle(paste0(cell, ': ROC Curve ', '(AUC = ', auc, ')')) +
  theme_minimal()
dev.off()
