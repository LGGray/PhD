library(CSCORE)
library(Seurat)
library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)

# Load the data
pbmc <- readRDS(commandArgs(trailingOnly = TRUE)[1])
DefaultAssay(pbmc) <- 'RNA'

cell <- levels(pbmc)[as.numeric(commandArgs(trailingOnly = TRUE)[2])]
# Subset to B cells from control and disease patients
pbmc_subset = pbmc[,pbmc$cellTypist %in% cell]
pbmc_subset_control <- pbmc_subset[, pbmc_subset$condition == "control"]
pbmc_subset_disease <- pbmc_subset[, pbmc_subset$condition == "disease"]

# Subsample the largest group to match the smallest group
if (nrow(pbmc_subset_control) > nrow(pbmc_subset_disease)) {
  pbmc_subset_control <- pbmc_subset_control[sample(1:nrow(pbmc_subset_control), nrow(pbmc_subset_disease)),]
} else {
  pbmc_subset_disease <- pbmc_subset_disease[sample(1:nrow(pbmc_subset_disease), nrow(pbmc_subset_control)),]
}

# Select top 5000 highly expressed genes
mean_exp = rowMeans(pbmc_subset@assays$RNA@counts/pbmc_subset$nCount_RNA)
genes_selected = names(sort.int(mean_exp, decreasing = T))[1:5000]

# Perform co-expression analysis on both conditions
CSCORE_control <- CSCORE(pbmc_subset_control, genes = genes_selected)
CSCORE_disease <- CSCORE(pbmc_subset_disease, genes = genes_selected)

# Extract co-expression estimates
coexp_control <- CSCORE_control$est
coexp_disease <- CSCORE_disease$est

# Obtain BH-adjusted p values
CSCORE_p_control <- CSCORE_control$p_value
p_matrix_BH = matrix(0, length(genes_selected), length(genes_selected))
p_matrix_BH[upper.tri(p_matrix_BH)] = p.adjust(CSCORE_p_control[upper.tri(CSCORE_p_control)], method = "BH")
p_matrix_BH_control <- p_matrix_BH + t(p_matrix_BH)
coexp_control[p_matrix_BH_control > 0.05] <- 0
save(coexp_control, file = paste0('cscore/', gsub("/|-| ", "_", cell), '_coexp_control.RData'))

CSCORE_p_disease <- CSCORE_disease$p_value
p_matrix_BH = matrix(0, length(genes_selected), length(genes_selected))
p_matrix_BH[upper.tri(p_matrix_BH)] = p.adjust(CSCORE_p_disease[upper.tri(CSCORE_p_disease)], method = "BH")
p_matrix_BH_disease <- p_matrix_BH + t(p_matrix_BH)
coexp_disease[p_matrix_BH_disease > 0.05] <- 0
save(coexp_disease, file = paste0('cscore/', gsub("/|-| ", "_", cell), '_coexp_disease.RData'))


# pdf('cscore/coexp_control.pdf')
# Heatmap(coexp_control, name = "Co-expression", 
# col = circlize::colorRamp2(c(-1, 0, 1), 
# c("blue", "white", "red")), 
# show_row_names = FALSE, show_column_names = FALSE)
# dev.off()

# pdf('cscore/coexp_disease.pdf')
# Heatmap(coexp_disease, name = "Co-expression",
# col = circlize::colorRamp2(c(-1, 0, 1),
# c("blue", "white", "red")),
# show_row_names = FALSE, show_column_names = FALSE)
# dev.off()


# Compute difference in co-expression between disease and control
coexp_diff <- coexp_disease - coexp_control

# Set up permutation test
set.seed(42)  # for reproducibility
n_permutations <- 100
permuted_diff <- matrix(0, nrow = length(genes_selected), ncol = n_permutations)

library(parallel)
cl <- makeCluster(8)
clusterExport(cl, varlist = c("pbmc_subset", "genes_selected", "CSCORE", "permuted_diff"))
permuted_diff <- parLapply(cl, 1:n_permutations, function(i) {
  # Randomly shuffle the group labels
  pbmc_perm <- pbmc_subset
  pbmc_perm$condition <- sample(pbmc_perm$condition)
  
  # Subset permuted data
  pbmc_perm_control <- pbmc_perm[, pbmc_perm$condition == "control"]
  pbmc_perm_disease <- pbmc_perm[, pbmc_perm$condition == "disease"]
  
  # Calculate co-expression for permuted data
  CSCORE_perm_control <- CSCORE(pbmc_perm_control, genes = genes_selected)
  CSCORE_perm_disease <- CSCORE(pbmc_perm_disease, genes = genes_selected)
  
  # Compute difference in permuted co-expression
  permuted_diff[, i] <- dim(CSCORE_perm_disease$est) - dim(CSCORE_perm_control$est)
})
stopCluster(cl)

save(coexp_diff, permuted_diff, file = "cscore/coexp_diff.RData")

# # Calculate p-values based on the null distribution of permuted differences
# p_values <- rowMeans(abs(coexp_diff) <= abs(permuted_diff))

# # Apply BH adjustment
# p_adjusted <- p.adjust(p_values, method = "BH")

# # Set co-expression entries with adjusted p-values > 0.05 to 0
# coexp_diff[p_adjusted > 0.05] <- 0

# # Continue with the WGCNA analysis using the significant differentially co-expressed pairs
# adj = WGCNA::adjacency.fromSimilarity(abs(coexp_diff), power = 1)
# TOM = WGCNA::TOMsimilarity(adj)
# dissTOM = 1 - TOM
# rownames(dissTOM) <- colnames(dissTOM) <- genes_selected
# hclust_dist = hclust(as.dist(dissTOM), method = "average")

# # Dynamic tree cutting for module identification
# memb = dynamicTreeCut::cutreeDynamic(dendro = hclust_dist, distM = dissTOM, deepSplit = 2, 
#                                      pamRespectsDendro = FALSE, minClusterSize = 10)
# names(memb) = genes_selected
# memb_tab <- table(memb)
# module_list = lapply(sort(unique(memb)), function(i_k) names(which(memb == i_k)))

