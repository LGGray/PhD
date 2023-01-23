library(boot)
library(dplyr)
library(tibble)
library(MASS)
library(caret)
library(pROC)
library(ggplot2)
source('~/R_code/functions/edgeR.list.R')
load('~/datasets/XCI/escapees.Rdata')

setwd('~/datasets/integrated/feature.selection')

args = commandArgs(trailingOnly=TRUE)
args <- as.numeric(args)
cell <- gsub('.HVG.X.RData', '', list.files('expr', pattern = '.HVG.X.RData'))[args]

load(paste0('expr/', cell,'.HVG.X.RData'))

# Set up k-fold cross validation
cv <- trainControl(method='cv', number=10)

# Split data into training and test sets 70:30
set.seed(42)
train_index <- createDataPartition(y=expr$condition, p=0.7, list=FALSE)
training_set <- expr[train_index, ]
testing_set <- expr[-train_index, ]

#fit a logistic regression model and use k-fold CV to evaluate performance
model <- train(condition ~., data = training_set, method = "glm", family = "binomial", trControl = cv)

importance <- varImp(model)$importance
features <- rownames(importance)[importance$Overall > mean(importance$Overall)]

# Remove low importance features and rebuild model
training_set <- training_set[,colnames(training_set) %in% c(features, 'condition')]
testing_set <- testing_set[,colnames(testing_set) %in% c(features, 'condition')]

model <- train(condition ~., data = training_set, method = "glm", family = "binomial", trControl = cv)

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


# plot ROC curve
predicted <- predict(model, newdata=testing_set, type='response')
rocobj <- roc(testing$condition, predicted)
auc <- round(auc(testing$condition, predicted),4)
pdf(paste0('feature.selection/logit/HVG.X', sub(' ', '_', cell), '.ROC.pdf'))
ggroc(rocobj, size=2) +
  ggtitle(paste0(cell, ': ROC Curve ', '(AUC = ', auc, ')')) +
  theme_minimal()
dev.off()

model <- list(logit.model=step.model, auc=auc)
save(model, file=paste0('logit/HVG.X/', sub(' ', '_', cell), '.model.RData'))
