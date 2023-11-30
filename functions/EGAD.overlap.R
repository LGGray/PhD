library(ggplot2)
library(reshape2)

# Overlap of EGAD results across studies
pSS <- read.delim('pSS_GSE157278/EGAD/GOBP/disease.pathway.txt')
UC <- read.delim('UC_GSE125527/EGAD/GOBP/disease.pathway.txt')
CD_colon <- read.delim('CD_Kong/colon/EGAD/GOBP/disease.pathway.txt')
CD_TI <- read.delim('CD_Kong/TI/EGAD/GOBP/disease.pathway.txt')
SLE <- read.delim('lupus_Chun/EGAD/GOBP/disease.pathway.txt')

celltypes <- lapply(list(pSS, UC, CD_colon, CD_TI, SLE), function(x) unique(x$celltype))
celltypes <- Reduce(intersect, celltypes)

studies <- list('pSS'=pSS, 'UC'=UC, 'CD_colon'=CD_colon, 'CD_TI'=CD_TI, 'SLE'=SLE)

# Measure overlap between test
jaccard_index <- function(list1, list2) {
  length(intersect(list1, list2)) / length(union(list1, list2))
}

# Create list to store jaccard scores
jaccard_lst <- list()
for(cell in celltypes){
    lst <- lapply(studies, function(x) subset(x, celltype == cell)$pathway)
    jaccard_scores <- combn(names(lst), 2, function(pair) {
        jaccard_index(lst[[pair[1]]], lst[[pair[2]]])
    })
    jaccard_lst[[cell]] <- jaccard_scores
}

result_lst <- list()
for(cell in celltypes){
    jaccard_matrix <- matrix(NA, nrow = length(names(studies)), ncol = length(names(studies)), 
    dimnames = list(names(studies), names(studies)))

    # Fill the upper triangle of the matrix with the Jaccard scores
    combinations <- combn(study_names, 2)
    for (i in 1:ncol(combinations)) {
        pair <- combinations[, i]
        jaccard_matrix[pair[1], pair[2]] <- jaccard_lst[[cell]][i]
    }
    result_lst[[cell]] <- jaccard_matrix
}

plot.data <- dplyr::bind_rows(lapply(result_lst, function(x) melt(x, na.rm = TRUE)), .id = 'celltype')
pdf('EGAD.celltype.jaccard.pdf')
ggplot(plot.data, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "blue", high = "red", na.value = "grey", limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(fill = "Jaccard Index", x = "", y = "") +
    facet_wrap(~celltype)
dev.off()

# Optionally, fill the lower triangle of the matrix (since it's symmetric)
jaccard_matrix[lower.tri(jaccard_matrix)] <- t(jaccard_matrix)[lower.tri(jaccard_matrix)]

# View the matrix
jaccard_matrix

# Converting to a matrix
jaccard_matrix <- matrix(jaccard_scores, nrow = length(studies), ncol = length(studies))
rownames(jaccard_matrix) <- names(studies)
colnames(jaccard_matrix) <- names(studies)
