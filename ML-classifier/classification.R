library(Seurat)
library(dplyr)
library(tibble)
library(e1071)
library(MASS)
library(nnet)
library(caret)
source('~/R_code/functions/edgeR.list.R')
load('~/datasets/XCI/escapees.Rdata')

setwd('datasets/integrated/')
pbmc <- readRDS('immune.combined.RDS')
edgeR <- edgeR.list('edgeR-LRT/', logfc=0.5)

cell = 'Treg'

xcape <- subset(edgeR[[cell]], gene %in% rownames(escape))
features <- xcape[order(xcape$FDR),][,1]

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

# Fit a logistical regression model to predict condition from all genes
glm.fits <- glm(
  condition ~ .,
  data = expr, family = binomial, subset = training
)
glm.fits
# Predicting condition from gene expression
glm.probs <- predict(glm.fits, testing,
                     type = "response")

glm.pred <- rep("Control", nrow(testing))
glm.pred[glm.probs > .5] = "Disease"

model_glm_pred <- ifelse(predict(glm.fits, type='response') > 0.5, "Yes", 'No')

calc_class_err <- function(actual, predicted){
  mean(actual != predicted)
}

calc_class_err(actual = training.data$condition, predicted = model_glm_pred)

table(predicted = model_glm_pred, actual = training.data$condition)

# Confusion matrix
table(glm.pred, testing$condition)
mean(glm.pred == testing$condition)

# Quadratric discriminant analysis - Assume gausian and each class has its own covariance matrix
qda.fit <- qda(condition ~ ., data = expr,
               subset = training)
# Check results
qda.fit

# Predict using testing data
qda.class <- predict(qda.fit, testing)$class
table(qda.class, testing$condition)
mean(qda.class == testing$condition)

# Naive Bayes
nb.fit <- naiveBayes(condition ~ ., data = expr,
                     subset = training)
# Check results
nb.fit

# Predict using testing data
nb.class <- predict(nb.fit, testing)
table(nb.class, testing$condition)
mean(nb.class == testing$condition)

## Multinomial regression
model_multi <- multinom(condition ~., data = training.data, trace=F)
summary(model_multi)
multi.pred <- predict(model_multi, newdata=testing)
table(multi.pred, testing$condition)
mean(multi.pred == testing$condition)

# testing which features have positive effect on model
features.new <- names(which(summary(model_multi)$coefficients > 0.1))
expr <- t(as.matrix(GetAssayData(pbmc.subset, assay='SCT')))
expr <- data.frame(expr[,which(colnames(expr) %in% features.new)])
expr$condition <- factor(pbmc.subset$condition, levels = c('Control', 'Disease'))

# Split data into training and test sets 80:20
training.data <- expr %>%
  rownames_to_column('cell') %>%
  group_by(condition) %>%
  sample_frac(0.8) %>%
  column_to_rownames('cell')
training <- rownames(expr) %in% rownames(training.data)

testing <- expr[!training,]

model_multi <- multinom(condition ~., data = training.data, trace=F)
summary(model_multi)
multi.pred <- predict(model_multi, newdata=testing)
table(multi.pred, testing$condition)
mean(multi.pred == testing$condition)

#######
full.model <- glm.fits <- glm(
  condition ~ .,
  data = expr, family = binomial, subset = training
)
coef(full.model)
step.model <- full.model %>% stepAIC(trace = F)
coef(step.model)

# Prediction accuracy of full model
probabilities <- full.model %>% predict(testing, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "Disease", "Control")
mean(predicted.classes == testing$condition)

# Prediction accuacy of stepwise model
probabilities <- predict(step.model, testing, type='response')
predicted.classes <- ifelse(probabilities > 0.5, "Disease", "Control")
mean(predicted.classes == testing$condition)

subset(xcape, gene %in% names(coef(step.model)))

