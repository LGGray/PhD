library(Seurat)
library(pvclust)

pbmc <- readRDS('onek1k.RDS')

# Determine sex of individuals by psuedobulked expression of chrY and XIST
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrY.Rdata')

PAR_genes <- c("PLCXD1", "GTPBP6", "LINC00685", "PPP2R3B", "FABP5P13", "KRT18P53", "SHOX", "RPL14P5", "CRLF2", "CSF2RA", "MIR3690", "RNA5SP498", "IL3RA", "SLC25A6", "LINC00106", "ASMTL-AS1", "ASMTL", "P2RY8", "AKAP17A", "ASMT", "DHRSX", "DHRSX-IT1", "ZBED1", "MIR6089", "CD99P1", "LINC00102", "CD99", "SPRY3", "DPH3P2", "VAMP7", "TRPC6P", "IL9R", "WASIR1", "WASH6P", "AJ271736.1", "DDX11L16")
chrY.nonPar <- chrY[!rownames(chrY) %in% PAR_genes,]
# pseduobulk expression matrix
exp <- AverageExpression(pbmc, assays='SCT', features=c('XIST', rownames(chrY)), group.by='individual')$SCT
exp <- AverageExpression(pbmc, assays='SCT', features=c('XIST', chrY.nonPar), group.by='individual')$SCT
exp <- AverageExpression(pbmc, assays='SCT', features=c('XIST', 'RPS4Y1'), group.by='individual')$SCT
exp <- AverageExpression(pbmc, assays='RNA', slot='counts', features=c('XIST', chrY.nonPar), group.by='individual')$RNA
exp <- scale(exp)

# First we infer sex based on expression of female specific XIST gene
XIST.expression <- exp[grep('XIST', rownames(exp)),]

# Perform hierarchical clustering to identify groups
dissimilarity <- dist(data.frame(t(exp)), method='euclidean')
cluster <- hclust(dissimilarity, method = 'centroid')

# Plot dendrogram
pdf('sex.dendrogram.pdf')
plot(cluster)
dev.off()

# K-means clustering on the hclust data
cluster.result <- cutree(cluster, k=2)
# Check the differencee between in expression of XIST between the two clusters
xist.1 <- XIST.expression[names(which(cluster.result==1))]
xist.2 <- XIST.expression[names(which(cluster.result==2))]
# Assign sex based on dendrogram
if(xist.1 > xist.2){
    sex.list <- ifelse(cluster.result == 1, 'F', 'M')
} else{
    sex.list <- ifelse(cluster.result == 1, 'M', 'F')
}
# Add sex to metadata
sle$sex.predict <- sex.list[pbmc$individual]

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

library(ggplot2)
# Plot the expression of XIST and RPS4Y1
ggplot(exp, aes(x='XIST', y='RPS4Y1')) + geom_point(aes(color=result$sex)) + theme_bw() + theme(legend.position='none')