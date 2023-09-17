# Determine sex of individuals by psuedobulked expression of chrY and XIST
predict.sex <- function(seurat.object, assay='RNA', slot='data', individual='individual'){
    
    # pseduobulk expression matrix
    exp <- AverageExpression(seurat.object, assays=assay, slot=slot, features=c('XIST', 'RPS4Y1'), group.by=individual)[[1]]
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
    xist.1 <- mean(XIST.expression[names(which(cluster.result==1))])
    xist.2 <- mean(XIST.expression[names(which(cluster.result==2))])
    # Assign sex based on dendrogram
    if(xist.1 > xist.2){
        sex.list <- ifelse(cluster.result == 1, 'F', 'M')
    } else{
        sex.list <- ifelse(cluster.result == 1, 'M', 'F')
    }
    # Add sex to metadata
    seurat.object$sex <- sex.list[pbmc$individual]
    return(seurat.object)
}

pbmc <- predict.sex(pbmc, assay='decontXcounts', slot='data', individual='individual')