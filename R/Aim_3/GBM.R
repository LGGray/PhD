library(Seurat)
library(dplyr)
library(tibble)
library(caret)
library(gbm)
library(pROC)
library(ggplot2)

setwd('~/datasets/integrated/feature.selection')

args = commandArgs(trailingOnly=TRUE)
args <- as.numeric(args)
cell <- gsub('.HVG.X.RData', '', list.files('expr', pattern = '.HVG.X.RData'))[args]

load(paste0('expr/', cell,'.HVG.X.RData'))

# Change condition to 0 1 
expr$condition <- ifelse(expr$condition == 'Control', 0, 1)

# Split data into training and test sets 70:30
set.seed(42)
train_index <- createDataPartition(y=expr$condition, p=0.7, list=FALSE)
training_set <- expr[train_index, ]
testing_set <- expr[-train_index, ]

# Method taken from http://uc-r.github.io/gbm_regression
hyper_grid <- expand.grid(
  shrinkage = c(.01, .1, .3),
  interaction.depth = c(1, 3, 5),
  n.minobsinnode = c(5, 10, 15),
  bag.fraction = c(.65, .8, 1), 
  optimal_trees = 0,              
  min_RMSE = 0                     
)

# grid search 
for(i in 1:nrow(hyper_grid)) {
  
  # reproducibility
  set.seed(123)
  
  # train model
  gbm.tune <- gbm(
    formula = condition ~ .,
    data = training_set,
    distribution = 'bernoulli',
    n.trees = 5000,
    interaction.depth = hyper_grid$interaction.depth[i],
    shrinkage = hyper_grid$shrinkage[i],
    n.minobsinnode = hyper_grid$n.minobsinnode[i],
    bag.fraction = hyper_grid$bag.fraction[i],
    train.fraction = .75,
    n.cores = 4, # will use all cores by default
    verbose = FALSE
  )
  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(gbm.tune$valid.error)
  hyper_grid$min_RMSE[i] <- sqrt(min(gbm.tune$valid.error))
}

parameters <- hyper_grid %>% 
  dplyr::arrange(min_RMSE) %>%
  head(1)

# train GBM model
gbm.fit.final <- gbm(
  formula = condition ~ .,
  data = training_set,
  distribution = 'bernoulli',
  shrinkage = parameters[1],
  interaction.depth = parameters[2],
  n.minobsinnode = parameters[3],
  bag.fraction = parameters[4], 
  n.trees = 5000,
  cv.folds = 10,
  n.cores = 4, # will use all cores by default
  verbose = FALSE
) 


# Determine importance of predictors and select best predictors
var.imp <- summary(gbm.fit, cBars=10, method=relative.influence, las=2)
features <- var.imp[var.imp[,2] > mean(var.imp[,2]),][,1]

# Remove low importance predictors
training_set <- training_set[,colnames(training_set) %in% c(features, 'condition')]
testing_set <- testing_set[,colnames(testing_set) %in% c(features, 'condition')]

# train GBM model
gbm.fit.final <- gbm(
  formula = condition ~ .,
  data = training_set,
  distribution = 'bernoulli',
  shrinkage = parameters[1],
  interaction.depth = parameters[2],
  n.minobsinnode = parameters[3],
  bag.fraction = parameters[4], 
  n.trees = 5000,
  cv.folds = 10,
  n.cores = 4, # will use all cores by default
  verbose = FALSE
)

# Determine the best number of trees
# best.iter <- gbm.perf(gbm.fit, method="cv")

# Refit model with best number of trees 
# fitControl = trainControl(method="cv", number=10, returnResamp = "all")
# gbm.fit2 = train(factor(condition)~., 
#                  data=training_set, 
#                  method="gbm",
#                  distribution="bernoulli", 
#                  trControl=fitControl, 
#                  verbose=F, 
#                  tuneGrid=data.frame(.n.trees=best.iter,
#                                      .shrinkage=0.01,
#                                      .interaction.depth=4,
#                                      .n.minobsinnode=1))

## Generate predictions
y_hats <- predict(
  
  ## Random forest object
  object=gbm.fit.final, 
  
  ## Data to use for predictions; remove the Species
  newdata=testing_set[, -(ncol(testing_set))])

## Print the accuracy
accuracy <- mean(y_hats == testing_set$condition)*100
print(paste('Accuracy on testing data: ', round(accuracy, 2), '%',  sep=''))
print(confusionMatrix(factor(testing_set$condition), y_hats))


#AUROC
predicted <- apply(predict(gbm.fit.final, newdata=testing_set, type='prob'), 1, max)
rocobj <- roc(testing_set$condition, predicted)
auc <- round(auc(testing_set$condition, predicted),4)
pdf(paste0('gbm/HVG.X/', sub(' ', '_', cell), '.ROC.pdf'))
ggroc(rocobj, size=2) +
  ggtitle(paste0(cell, ': ROC Curve ', '(AUC = ', auc, ')')) +
  theme_minimal()
dev.off()

model <- list(gbm.model=gbm.fit2, auc=auc)
save(model, file=paste0('gbm/HVG.X/', sub(' ', '_', cell), '.model.RData'))


