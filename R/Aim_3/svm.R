library(dplyr)
library(tibble)
library(caret)
library(e1071)
library(pROC)
library(ggplot2)
library(ROCR)

setwd('~/datasets/integrated/feature.selection')

args = commandArgs(trailingOnly=TRUE)
args <- as.numeric(args)
cell <- gsub('.deg.RData', '', list.files('expr', pattern = '.deg.RData'))[args]

load(paste0('expr/', cell,'.deg.RData'))

# Split data into training and test sets 70:30
set.seed(42)
train_index <- createDataPartition(y=expr$condition, p=0.7, list=FALSE)
training_set <- expr[train_index, ]
testing_set <- expr[-train_index, ]

# train svm 
tune.out <- tune(svm, condition ~ ., data = training_set,
                 kernel = 'radial',
                 ranges = list(
                   cost = c(0.1, 1, 10, 100, 1000),
                   gamma = c(0.5, 1, 2, 3, 4)
                 ))
print(summary(tune.out))

svm.model <- tune.out$best.model

print(table(
  true = testing_set$condition,
  pred = predict(
    svm.model, newdata = testing_set
  )
))

rocplot <- function(pred, truth, ...){
  predob <- prediction(pred, truth)
  perf <- performance(predob, 'tpr', 'fpr')
                      plot(perf, ...)
}
fitted <- attributes(
  predict(svm.model, 
          testing_set, decision.values = T))$decision.values
pred <- prediction(fitted, testing_set$condition)
auc.perf = performance(pred, measure = "auc")
auc.perf@y.values[[1]]
pdf(paste0('svm/deg/', cell, '.AUROC.pdf'))
rocplot(fitted, testing_set$condition, col='red', AUC=T)
mtext(paste0(cell, ': ROC Curve (AUC = ', round(auc.perf@y.values[[1]], 2), ')'))
dev.off()

model <- list(svm.model=svm.model, auc=round(auc.perf@y.values[[1]], 2))
save(model, file=paste0('svm/deg/', cell, '.model.RData'))
