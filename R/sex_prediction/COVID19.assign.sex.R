library(Seurat)
library(transformGamPoi)

# Setwd 
setwd('/directflow/SCCGGroupShare/projects/lacgra/seurat.object/COVID19')

# Read in Seurat object
pbmc <- readRDS('local.rds')

print(pbmc)

# Normalise data with Delta method-based variance stabilizing
exp.matrix <- GetAssayData(pbmc, slot = 'counts')
exp.matrix.transformed <- acosh_transform(exp.matrix)

# Add transformed data to Seurat object
pbmc <- SetAssayData(object=pbmc, slot = 'counts', new.data=exp.matrix.transformed)



# Read in Y chromosome genes
# load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrY.Rdata')
# Remove PAR genes from chrY
# PAR_genes <- c("PLCXD1", "GTPBP6", "LINC00685", "PPP2R3B", "FABP5P13", "KRT18P53", "SHOX", "RPL14P5", "CRLF2", "CSF2RA", "MIR3690", "RNA5SP498", "IL3RA", "SLC25A6", "LINC00106", "ASMTL-AS1", "ASMTL", "P2RY8", "AKAP17A", "ASMT", "DHRSX", "DHRSX-IT1", "ZBED1", "MIR6089", "CD99P1", "LINC00102", "CD99", "SPRY3", "DPH3P2", "VAMP7", "TRPC6P", "IL9R", "WASIR1", "WASH6P", "AJ271736.1", "DDX11L16")
# chrY.nonPar <- chrY[!rownames(chrY) %in% PAR_genes,]
# Determine sex of individuals by psuedobulked expression of chrY.nonPar and XIST
exp <- AverageExpression(pbmc, assays='RNA', slot='counts', features=c('XIST', 'RPS4Y1'), group.by='individual')$RNA
exp <- scale(exp)

save(exp, file='sexpredict.exp.Rdata')

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
