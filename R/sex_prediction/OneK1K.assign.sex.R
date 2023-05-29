library(Seurat)

setwd('/directflow/SCCGGroupShare/projects/lacgra/seurat.object')

# Read in Seurat object
pbmc <- readRDS('onek1k.RDS')

# Determine sex of individuals by psuedobulked expression of chrY.nonPar and XIST
exp <- AverageExpression(pbmc, assays='SCT', slot='counts', features=c('XIST', 'RPS4Y1'), group.by='individual')$SCT
exp <- scale(exp)

# First we infer sex based on expression of female specific XIST gene as sanity check
XIST.expression <- exp[grep('XIST', rownames(exp)),]

# Perform hierarchical clustering to identify groups
set.seed(42)
dissimilarity <- dist(data.frame(t(exp)), method='euclidean')
cluster <- hclust(dissimilarity, method = 'median')

save(dissimilarity, file='sexpredict.dissimilarity.Rdata')
save(cluster, file='sexpredict.cluster.Rdata')

# Plot dendrogram
pdf('sex.dendrogram.pdf')
plot(cluster, labels=FALSE, main = 'SLE Dendrogram (M | F)', sub = 'RPS4Y1 | XIST')
dev.off()

# K-means clustering on the hclust data
cluster.result <- cutree(cluster, k=2)
# Check the differencee between in expression of XIST between the two clusters
xist.1 <- mean(XIST.expression[names(which(cluster.result==1))])
xist.2 <- mean(XIST.expression[names(which(cluster.result==2))])
# Assign sex based on dendrogram
if(xist.1 > xist.2){
    sex.list <- ifelse(cluster.result == 1, 'F', 'M')
} else{
    sex.list <- ifelse(cluster.result == 1, 'M', 'F')
}
# Add sex to metadata
pbmc$sex.predict <- sex.list[pbmc$individual]

result <- unique(pbmc@meta.data[,c('individual', 'sex', 'sex.predict')])

# Calculate the number of correct predictions
correct_predictions <- sum(result$sex == result$sex.predict)
# Calculate the accuracy
accuracy <- correct_predictions / nrow(result)
# Print the accuracy
cat("Accuracy:", accuracy)

# Create a table of the predicted and known sex
confusion_matrix <- table(result$sex.predict, result$sex)
# Calculate the number of true positive (TP), false positive (FP), true negative (TN), and false negative (FN) predictions
TP <- confusion_matrix[2, 2]
FP <- confusion_matrix[2, 1]
TN <- confusion_matrix[1, 1]
FN <- confusion_matrix[1, 2]
# Calculate the sensitivity (recall)
sensitivity <- TP / (TP + FN)
# Calculate the specificity
specificity <- TN / (TN + FP)
# Calculate the precision
precision <- TP / (TP + FP)
# Calculate the F1-score
F1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)
# Print the confusion matrix and the performance metrics
cat("Confusion matrix:\n")
print(confusion_matrix)
cat("\nSensitivity:", sensitivity)
cat("\nSpecificity:", specificity)
cat("\nPrecision:", precision)
cat("\nF1-score:", F1_score)

# Create a heatmap of the expression of XIST and chrY.nonPar genes
library(ggplot2)
library(reshape2)

exp.melt <- melt(exp)
pdf('sex.predict.heatmap.pdf')
ggplot(exp.melt, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + 
    scale_fill_gradient(low="white", high="red") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    theme(axis.text.y = element_blank()) +
    labs(x='Individual',y='', fill='z-scored Expression')
dev.off()

# Load the PRROC package
library(PRROC)
# Extract the expression values of XIST and RPS4Y1 genes
XIST <- exp["XIST", ]
RPS4Y1 <- exp["RPS4Y1", ]
# Convert the known sex to a binary factor
sex <- ifelse(result$sex == "F", 1, 0)
# Calculate the prROC curve for XIST
prroc.XIST <- pr.curve(sex, XIST, curve=T)
# Calculate the prROC curve for RPS4Y1
prroc.RPS4Y1 <- pr.curve(sex, RPS4Y1, curve=T)
# Plot the ROC curves
pdf('sex.predict.prROC.pdf')
plot(prroc.XIST, col = "blue", main = "ROC curves for XIST and RPS4Y1")
plot(prroc.RPS4Y1, col = "red", add = TRUE)
legend("bottomright", c("XIST", "RPS4Y1"), col = c("blue", "red"), lty = 1)
dev.off()