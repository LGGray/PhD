library(ggplot2)
library(dplyr)
library(Seurat)
library(car)
library(rstatix)
library(glmGamPoi)
library(reshape2)
library(lmerTest)
library(emmeans)

if(dir.exists('differential.variance') != TRUE){dir.create('differential.variance')}

pbmc <- readRDS('pbmc.female.RDS')

result <- lapply(levels(pbmc), function(cell){
    # Select cell type
    print(cell)
    # subset object by cell type
    pbmc.cell <- subset(pbmc, cellTypist == cell)
    # check if there are enough cell in both conditions and skip if not
    if(length(unique(pbmc.cell$condition)) != 2){
        return("Not enough conditions")
        next
    }
    # Keep genes with expression in 5% of cells
    keep <- rowSums(pbmc.cell@assays$RNA@counts > 0) > ncol(pbmc.cell) * 0.05
    features <- names(keep[keep == T])
    pbmc.cell <- subset(pbmc.cell, features=features)

    individuals <- unique(pbmc.cell$individual)
    res <- data.frame(matrix(NA, nrow=nrow(pbmc.cell), ncol=length(individuals)))
    for(i in 1:length(individuals)){
        counts <- GetAssayData(subset(pbmc.cell, individual==individuals[i]), slot='counts')
        fit <- glm_gp(as.matrix(counts), size_factors = FALSE, verbose = F)
        res[,i] <- fit$overdispersions
    }
    rownames(res) <- rownames(pbmc.cell)
    colnames(res) <- individuals
    res <- cbind(gene=rownames(res), res)

    # Test for differential variance
    # Melt the data
    res.melt <- melt(res)
    # Rename the variables
    names(res.melt) <- c("gene", "individual", "value")
    # Create a group variable based on the column names
    res.melt$group <- ifelse(grepl("HC", res.melt$individual), "control", "disease")

    res.melt.chrX <- subset(res.melt, gene %in% rownames(chrX))

    model <- lmer(value ~ group + (1 | individual), data = res.melt.chrX)
    summary(model)
    
    emmeans(model, pairwise ~ group)

    # Perform wilcox test for each gene
    test.result <- res.melt %>%
    group_by(gene) %>%
    summarise(
        statistic = wilcox.test(value ~ group)$statistic,
        p.value = wilcox.test(value ~ group)$p.value
    ) %>%
    mutate(FDR = p.adjust(p.value, method = "fdr")) %>%
    as.data.frame()

    # summarize the results
    print(summary(test.result$FDR))

    subset(test.result, FDR < 0.05)

    write.table(test.result, file=paste0('differential.variance/', gsub(' |/|-', '_', cell), '.txt'), sep='\t', row.names=F, quote=F)
})


