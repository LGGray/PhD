library(dplyr)
library(tibble)
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
#source('~/R_code/functions/edgeR.list.R')
# load('~/datasets/XCI/escapees.Rdata')
# load('~/datasets/XCI/chrX.Rdata')

setwd('~/datasets/integrated/feature.selection')

args = commandArgs(trailingOnly=TRUE)
args <- as.numeric(args)
cell <- gsub('.deg.RData', '', list.files('expr', pattern = 'deg.RData'))[args]

load(paste0('expr/', cell,'.deg.RData'))

# Set up k-fold cross validation
cv <- trainControl(method='cv', number=10)

# Split data into training and test sets 70:30
set.seed(42)
train_index <- createDataPartition(y=expr$condition, p=0.7, list=FALSE)
training_set <- expr[train_index, ]
testing_set <- expr[-train_index, ]

## Train a random forest model
forest <- train(
  # Formula. We are using all variables to predict condition
  condition~., 
  # Source of data;
  data=training_set, 
  # `rf` method for random forest
  method='rf',
  # Add repeated cross validation as trControl
  trControl=cv,
  # Accuracy to measure the performance of the model
  metric='Accuracy')

## Print out the details about the model
print(forest$finalModel)

## Get variable importance, and turn into a data frame
var_imp <- varImp(forest, scale=FALSE)$importance
var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)
## Select features with importance > mean(importance)
new.features <- var_imp[var_imp$importance > mean(var_imp$importance),]$variables

# Expr with new features
expr <- data.frame(expr[,colnames(expr) %in% new.features], condition = expr$condition)

# Split data into training and test sets 70:30
set.seed(42)
train_index <- createDataPartition(y=expr$condition, p=0.7, list=FALSE)
training_set <- expr[train_index, ]
testing_set <- expr[-train_index, ]

## Train a random forest model
forest <- train(
  # Formula. We are using all variables to predict condition
  condition~., 
  # Source of data;
  data=training_set, 
  # `rf` method for random forest
  method='rf',
  # Add repeated cross validation as trControl
  trControl=cv,
  # Accuracy to measure the performance of the model
  metric='Accuracy')

## Print out the details about the model
print(forest$finalModel)

## Generate predictions
y_hats <- predict(
  
  ## Random forest object
  object=forest, 
  
  ## Data to use for predictions; remove the Species
  newdata=testing_set[,-ncol(testing_set)])

## Print the accuracy
accuracy <- mean(y_hats == testing_set$condition)*100
print(paste('Accuracy on testing data: ', round(accuracy, 2), '%',  sep=''))
print(confusionMatrix(testing_set$condition, y_hats))

# plot ROC curve
predicted <- apply(predict(forest, newdata=testing_set, type='prob'), 1, max)
rocobj <- roc(testing_set$condition, predicted)
auc <- round(auc(testing_set$condition, predicted),4)
pdf(paste0('randomForest/deg/', cell, '.AUROC.pdf'))
ggroc(rocobj, size=2) +
  ggtitle(paste0(cell, ': ROC Curve ', '(AUC = ', auc, ')')) +
  theme_minimal()
dev.off()

model <- list(gbm.model=forest, auc=auc)
save(model, file=paste0('randomForest/deg/', cell, '.model.RData'))
