library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)

pbmc <- readRDS(commandArgs(trailingOnly=TRUE))
metadata <- pbmc@meta.data
auc_mtx <- read.csv('SCENIC/SCENIC.auc.csv', row.names=1)
colnames(auc_mtx) <- gsub('\\.', '', colnames(auc_mtx))
# Match cellIDs to pbmc
common_rows <- intersect(rownames(metadata), rownames(auc_mtx))
metadata <- metadata[rownames(metadata) %in% common_rows, ]
auc_mtx <- auc_mtx[rownames(auc_mtx) %in% common_rows, ]

# Add condition and celltype to auc_mtx
auc_mtx$condition <- metadata$condition[match(rownames(metadata), rownames(auc_mtx))]
auc_mtx$celltype <- metadata$cellTypist[match(rownames(metadata), rownames(auc_mtx))]
auc_mtx$individual <- metadata$individual[match(rownames(metadata), rownames(auc_mtx))]

# reshape data for plotting and statistical analysis
auc_mtx_melt <- reshape2::melt(auc_mtx, id.vars = c('condition', 'celltype', 'individual'))
auc_mtx_melt$condition <- factor(auc_mtx_melt$condition, levels = c('disease', 'control'))

# Create TF_celltype column
auc_mtx_melt$TF_celltype <- paste(auc_mtx_melt$variable, auc_mtx_melt$celltype, sep = ':')

test <- split(auc_mtx_melt, auc_mtx_melt$TF_celltype)[[1]]
wilcox.test(test[test$condition == 'disease', 'value'], test[test$condition == 'control', 'value'])


auc_sum <- aggregate(auc_mtx_melt$value, by = list(auc_mtx_melt$condition, auc_mtx_melt$individual, auc_mtx_melt$celltype, auc_mtx_melt$variable), FUN = sum)
colnames(auc_sum) <- c("condition", "individual", "celltype", "TF", "AUC")
auc_sum$celltype_TF <- paste(auc_sum$celltype, auc_sum$TF, sep = '_')

auc_avg <- aggregate(auc_mtx_melt$value, by = list(auc_mtx_melt$condition, auc_mtx_melt$individual, auc_mtx_melt$celltype, auc_mtx_melt$variable), FUN = mean)
colnames(auc_avg) <- c("condition", "individual", "celltype", "TF", "AUC")
auc_avg$celltype_TF <- paste(auc_avg$celltype, auc_avg$TF, sep = '_')

# Split by celltype_TF
auc_sum_split <- split(auc_sum, auc_sum$celltype_TF)
auc_avg_split <- split(auc_avg, auc_avg$celltype_TF)

# load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')

# reg <- read.csv('SCENIC/reg.csv', skip=3, header=F)
# reg_greater <- subset(reg, V1 %in% unique(wilcox_greater$TF))

# gene_targets <- lapply(reg_greater$V9, function(x) gsub("'", '', unlist(stringr::str_extract_all(x, "'[^']+'"))))
# names(gene_targets) <- reg_greater$V1
# gene_targets_escape <- lapply(gene_targets, function(x) x[x %in% rownames(escape)])
# gene_targets_escape <- gene_targets_escape[sapply(gene_targets_escape, length) > 0]

set.seed(123)
results <- list()
for(x in 1:length(auc_sum_split)){
  control_values <- auc_sum_split[[x]][auc_sum_split[[x]]$condition == 'control','AUC']
  disease_values <- auc_sum_split[[x]][auc_sum_split[[x]]$condition == 'disease', 'AUC']

  # Calculate actual difference in means
  actual_diff <- mean(disease_values) - mean(control_values)

  # Initialize a null distribution
  null_distribution <- numeric(10000)

  # Perform permutations
  for (i in 1:10000) {
      # Shuffle the AUC values
      shuffled_values <- sample(c(control_values, disease_values))
      
      # Split into new groups
      perm_control <- shuffled_values[1:length(control_values)]
      perm_disease <- shuffled_values[(length(control_values)+1):length(shuffled_values)]
      
      # Calculate the difference in means for the permuted groups
      null_distribution[i] <- mean(perm_disease) - mean(perm_control)
  }

  # Calculate the p-value
  p_value <- sum(null_distribution >= actual_diff) / length(null_distribution)
  results[[x]] <- c(actual_diff, p_value)
}

results <- data.frame(do.call(rbind, results))
colnames(results) <- c('actual_diff', 'p_value')
results <- cbind(celltype_TF = names(auc_sum_split), results)
results$celltype <- gsub(':.+', '', results$celltype_TF)
results$TF <- gsub('.+:', '', results$celltype_TF)

# results %>% 
# group_by(celltype) %>% 
# reframe(FDR = p.adjust(p_value, method = "fdr"))

results$FDR <- p.adjust(results$p_value, method = "fdr")
results <- results[,c('celltype_TF', 'celltype', 'TF', 'actual_diff', 'p_value', 'FDR')]
write.csv(results, 'SCENIC/permutation_test.csv', row.names = F)

### Null hypothesis testing for average AUC ###
set.seed(123)
results_avg <- list()
for(x in 1:length(auc_avg_split)){
  control_values <- auc_avg_split[[x]][auc_avg_split[[x]]$condition == 'control','AUC']
  disease_values <- auc_avg_split[[x]][auc_avg_split[[x]]$condition == 'disease', 'AUC']

  # Calculate actual difference in means
  actual_diff <- mean(disease_values) - mean(control_values)

  # Initialize a null distribution
  null_distribution <- numeric(10000)

  # Perform permutations
  for (i in 1:10000) {
      # Shuffle the AUC values
      shuffled_values <- sample(c(control_values, disease_values))
      
      # Split into new groups
      perm_control <- shuffled_values[1:length(control_values)]
      perm_disease <- shuffled_values[(length(control_values)+1):length(shuffled_values)]
      
      # Calculate the difference in means for the permuted groups
      null_distribution[i] <- mean(perm_disease) - mean(perm_control)
  }

  # Calculate the p-value
  p_value <- sum(null_distribution >= actual_diff) / length(null_distribution)
  results_avg[[x]] <- c(actual_diff, p_value)
}

results_avg <- data.frame(do.call(rbind, results_avg))
colnames(results_avg) <- c('actual_diff', 'p_value')
results_avg <- cbind(celltype_TF = names(auc_avg_split), results_avg)
results_avg$celltype <- gsub(':.+', '', results_avg$celltype_TF)
results_avg$TF <- gsub('.+:', '', results_avg$celltype_TF)

results_avg$FDR <- p.adjust(results_avg$p_value, method = "fdr")
results_avg <- results_avg[,c('celltype_TF', 'celltype', 'TF', 'actual_diff', 'p_value', 'FDR')]
write.csv(results_avg, 'SCENIC/permutation_test_avg.csv', row.names = F)

