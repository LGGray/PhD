###--------
# A package to determine biological sex from transcriptomic data 
# The aim of this package is to assist in assigning sex to bulk RNA sequencing 
# and single cell RNA sequencing datasets downloaded from online databases.
# PCA code from https://tavareshugo.github.io/data-carpentry-rnaseq/03_rnaseq_pca.html

library(ggplot2)
library(reshape2)
library(edgeR)
library(tidyverse)

load('~/external/ClusterHome/datasets/XCI/escapees.Rdata')
load('~/external/ClusterHome/datasets/XCI/chrX.Rdata')

bulkRNA <- function(infile){
  
}

scRNA <- function(seurat.object){
  
}

funcname <- function(bulk=TRUE, infile){
  if(bulk){bulkRNA(infile)
  }else{
      scRNA(infile)
    }
}

infile <- read.csv('~/Dropbox (Garvan)/LGray PhD/SciX 2022 files/HD.data.csv', row.names = 1)

# create metadata
metadata <- data.frame(sample=colnames(infile))
metadata$condition <- gsub('[0-9]', '', read.delim('~/Dropbox (Garvan)/LGray PhD/SciX 2022 files/HD.metadata.txt')$condition)
XIST.counts <- infile[grep('XIST', rownames(infile)),]
metadata$sex <- apply(XIST.counts, 1, function(x) ifelse(x> 20, 'F', 'M'))

# TMM normalisation
y = DGEList(counts=infile)
tmm <- cpm(y)

### PCA Clustering -----------------------------
# Remove empty rows for PCA
pca.data <- tmm.subset[rowSums(tmm.subset) > 0,]
res.pca <- prcomp(t(pca.data), scale = T)

ggplot(data.frame(res.pca$x), aes(x=PC1, y=PC2, colour=metadata$sex)) + geom_point() +
  geom_text(data=data.frame(res.pca$x), aes(x=PC1, y=PC2, label=rownames(res.pca$x)))
  
pc_eigenvalues <- res.pca$sdev^2
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)),
                             variance = pc_eigenvalues) %>%
                             mutate(pct = variance/sum(variance)*100) %>%
                             mutate(pct_cum = cumsum(pct))
# Scree Plot
pc_eigenvalues %>%
  ggplot(aes(x=PC)) +
  geom_col(aes(y=pct))+
  geom_line(aes(y=pct_cum, group = 1)) +
  geom_point(aes(y = pct_cum)) +
  labs(x = 'PC', y = 'Fraction variance explained')

# PCA plot
pc_scores <- res.pca$x
pc_scores <- pc_scores %>% as_tibble(rownames = "sample")
pc_scores %>% 
  ggplot(aes(x = PC1, y = PC2, colour=metadata$condition)) +
  geom_point(shape=metadata$sex, size=5)

### Clustering on escape genes ----------------------
tmm.subset <- subset(tmm, gsub('ENSG[0-9]+\\.[0-9]+\\|', '', rownames(tmm)) %in% rownames(chrY))
rownames(tmm.subset) <- gsub('ENSG[0-9]+\\.[0-9]+\\|', '', rownames(tmm.subset))

# Scale data 
hclust_matrix <- tmm.subset %>%
  t() %>%
  scale() %>%
  t()

# Calculate distance
gene_dist <- dist(hclust_matrix)

# Fill in NA with 0
gene_dist[is.na(gene_dist)] <- 0

# Heirarchical clustering
gene_hclust <- hclust(t(gene_dist), method = "complete")
plot(gene_hclust, labels = F)
abline(h=6, col='red', lwd=2)

# Cut dendogram into n groups
cutree(gene_hclust, k=2)

gene_cluster <- cutree(gene_hclust, k=2) %>%
  enframe() %>%
  rename(gene = name, cluster = value)
head(gene_cluster)





