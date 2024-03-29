library(ggplot2)
library(dplyr)
library(Seurat)
library(glmGamPoi)
library(reshape2)
library(qvalue)

if(dir.exists('differential.variance') != TRUE){dir.create('differential.variance')}

pbmc <- readRDS('pbmc.female.RDS')

lapply(levels(pbmc), function(cell){
    
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
    res.melt$group <- ifelse(grepl("HC", res.melt$individual), 0, 1)
    res.melt$group <- factor(res.melt$group)
    # Add batch effect
    batch <- tapply(pbmc.cell$SV1, pbmc.cell$individual, sum)
    res.melt$batch <- batch[res.melt$individual]

    # Fit a linear model to each gene
    result <- matrix(NA, nrow=nrow(res), ncol=2)
    for(i in 1:nrow(res)){
        model <- lm(value ~ batch + group, data = subset(res.melt, gene %in% res$gene[i]))
        p.value <- summary(model)$coefficients[3,4]
        result[i,1] <- p.value
        logFC <- coef(model)[3]
        result[i,2] <- logFC
    }
    result <- data.frame(logFC=result[,2], p.value=result[,1])
    result <- cbind(gene=res$gene, result)
    result$FDR <- qvalue(p = result$p.value)$qvalues

    # # Perform wilcox test for each gene
    # wilcox.result <- res.melt %>%
    # group_by(gene) %>%
    # summarise(
    #     statistic = wilcox.test(value ~ group, alternative='greater')$statistic,
    #     p.value = wilcox.test(value ~ group)$p.value
    # ) %>%
    # mutate(FDR = qvalue(p = result$p.value)$qvalues) %>%
    # as.data.frame()

    write.table(result, file=paste0('differential.variance/', gsub(' |/|-', '_', cell), '.txt'), sep='\t', row.names=F, quote=F)
})