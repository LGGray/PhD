library(ComplexHeatmap)
library(circlize)

load('new_pseudobulk/figures/selected_features.RData')

celltypes <- c('Memory.B.cells.chrX', 'Regulatory.T.cells.chrX', 'Non.classical.monocytes.chrX',
'Tcm.Naive.helper.T.cells.chrX', 'pDC.chrX', 'CD16+.NK.cells.chrX')

for(i in 1:length(celltypes)){
    exp <- readRDS(paste0('../', celltypes[i], '.RDS'))
    exp <- exp[,c('class', selected_features$all_features[[celltypes[i]]])]
    exp <- exp[order(exp$class),]
    class <- exp$class
    exp <- t(scale(as.matrix(exp[,-1])))
    dend <- cluster_within_group(exp, class)

    # Find optimal k for k-means clustering
    kmeans <- list()
    for (j in 2:10) {
    kmeans[[j]] <- kmeans(exp, centers = j)
    }
    wss <- sapply(kmeans, function(km) sum(km$withinss))
    # select lowest wss
    k <- which.min(wss)+1

    pdf(paste0('expression_heatmaps/', celltypes[i], '.pdf'))
    anno <- HeatmapAnnotation(df = data.frame(class = class), col = list(class = c('disease' = 'orange', 'control' = 'purple')))
    plot <- Heatmap(exp, col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")), row_names_gp = gpar(fontsize = 3), 
    name = 'z-score', show_row_names = TRUE, show_column_names = FALSE, column_title = celltypes[i],
    top_annotation = anno, cluster_rows = TRUE, cluster_columns = dend, row_km = 2)
    print(plot)
    dev.off()
}