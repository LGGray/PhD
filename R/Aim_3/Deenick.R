library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(dplyr)
library(UpSetR)

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')
# feature.files <- list.files('psuedobulk/ML.models/ensemble/features/', pattern='perm.*.chrX.txt', full.names=TRUE)
# feature.files <- feature.files[c(11,12,13,17,18,19,20,21,22)]
# features <- lapply(feature.files, read.delim, header=T)
# names(features) <- replace.names(gsub('perm.|.chrX.txt', '', basename(feature.files)))

top_models <- c('Tcm.Naive.helper.T.cells', 'Regulatory.T.cells', 'Non.classical.monocytes', 'Memory.B.cells')
feature.list <- lapply(top_models, function(x){
    read.csv(paste0('psuedobulk/features/combined_features.', x, '.chrX.csv'))[,1]
})
names(feature.list) <- replace.names(top_models)


nichenet <- readRDS('/directflow/SCCGGroupShare/projects/lacgra/NicheNet/lr_network_human.RDS')

surface.features <- lapply(feature.list, function(x){
    x[x %in% unique(c(nichenet$from, nichenet$to))]
})

surface.mtx <- fromList(surface.features)
rownames(surface.mtx) <- unique(unlist(surface.features))
pdf('psuedobulk/ML.plots/Selected.surface.features.heatmap.pdf')
Heatmap(as.matrix(surface.mtx), show_row_names=T, show_column_names=T, show_heatmap_legend = FALSE)
dev.off()

pbmc <- readRDS('pbmc.female.control-managed.RDS')

infile <- 'psuedobulk/Regulatory.T.cells.chrX.RDS'
cell <- replace.names(gsub('.chrX.RDS', '', basename(infile)))
expr <- readRDS(infile)
expr <- expr[order(expr$class),]
Condition <- expr$class
names(Condition) <- expr$individual
expr <- expr[,feature.list[[cell]]]
expr <- t(scale(expr))

# Heatmap of selected features
pdf(paste0('psuedobulk/ML.plots/', gsub('.chrX.RDS', '', basename(infile)), '.chrX.heatmap.pdf'))
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
column_ha = HeatmapAnnotation(Condition=Condition, col=list(Condition=c('control'='purple', 'disease'='pink')))
Heatmap(expr, name='z-score', column_title=cell, show_row_names=T, show_column_names=F, 
show_heatmap_legend = FALSE, top_annotation = column_ha, col=col_fun,
clustering_distance_columns = 'euclidean', clustering_method_columns = 'complete',
clustering_distance_rows = 'euclidean', clustering_method_rows = 'complete',
column_split = Condition)
dev.off()

# Boxplot of selected features
expr.melt <- reshape2::melt(expr, varnames = c("gene", "individual"), value.name = "zscore")
expr.melt$Condition <- factor(Condition[match(expr.melt$individual, names(Condition))])
pdf(paste0('psuedobulk/ML.plots/', gsub('.chrX.RDS', '', basename(infile)), '.chrX.boxplot.pdf'))
ggplot(expr.melt, aes(x=gene, y=zscore, fill=Condition)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),
    plot.margin = unit(c(1,1,1,3), "lines")) +
    labs(y='z-score', x='') + ggtitle(cell)
dev.off()

# Wilcoxon rank sum test
wilcox <- lapply(split(expr.melt, expr.melt$gene), function(x){
    wilcox.test(x$zscore ~ x$Condition)$p.value
})
wilcox.results <- data.frame(gene=names(wilcox), p.value=unlist(wilcox))
wilcox.results$FDR <- p.adjust(wilcox.results$p.value, method='BH')
wilcox.results[order(wilcox.results$FDR),]

# Correlation matrix of selected features
cor.results <- cor(t(expr), method='spearman')
# Mask lower triangle with NA
cor.results[lower.tri(cor.results)] <- 0
pdf(paste0('psuedobulk/ML.plots/', gsub('.chrX.RDS', '', basename(infile)), '.chrX.cor.pdf'))
Heatmap(cor.results, name='correlation', column_title=cell, show_row_names=T, show_column_names=T,
show_heatmap_legend = FALSE, col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
cluster_columns = FALSE, cluster_rows=FALSE)
dev.off()




meta <- data.frame(unique(pbmc@meta.data[,c('individual', 'condition')]))

# Plot pseudobulke expression of selected features across cell types
for(gene in surface.features){
    psuedobulk_list <- list()
    for(cell in levels(pbmc)){
        pbmc.subset <- subset(pbmc, idents=cell)
        psuedobulk <- data.frame(t(AggregateExpression(pbmc.subset, group.by='individual', features=gene)$RNA))
        colnames(psuedobulk) <- 'value'
        psuedobulk$condition <- meta[match(rownames(psuedobulk), meta$individual),'condition']
        rownames(psuedobulk) <- NULL
        psuedobulk$value <- scale(psuedobulk$value)
        psuedobulk_list[[cell]] <- psuedobulk
    }

    psuedobulk <- bind_rows(psuedobulk_list, .id='celltype')
    psuedobulk$condition <- factor(psuedobulk$condition, levels=c('control', 'disease'))

    write.table(psuedobulk, paste0('psuedobulk/target.celltype.expression/', gene, '.txt'), sep='\t', row.names=F, quote=F)

    quantiles <- quantile(psuedobulk$value, probs = c(0.01, 0.99), na.rm = TRUE)
    pdf(paste0('Deenick/',gene,'.boxplot.pdf'))
    print(ggplot(psuedobulk, aes(x=celltype, y=value, fill=condition)) +
        geom_boxplot(outlier.shape = NA) +
        theme_bw() +
        coord_cartesian(ylim = c(quantiles[1], quantiles[2])) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),
        plot.margin = unit(c(1,1,1,3), "lines")) +
        labs(y='z-score', x='') + ggtitle(gene))
    dev.off()
}

# Plot percentage of cells expressing selected features across cell types
for(gene in surface.features){
    gene_expression <- FetchData(pbmc, vars = gene)
    gene_expression$individual <- pbmc$individual
    gene_expression$condition <- pbmc$condition
    gene_expression$celltype <- pbmc$cellTypist
    gene_sym <- sym(gene)
    perc.expression <- gene_expression %>%
    group_by(individual, celltype) %>%
    summarise(
        NumCellsExpressing = sum(!!gene_sym > 0)/n(), 
        .groups = 'drop'
    ) %>%
    left_join(meta, by = "individual") %>%
    data.frame()

    pdf(paste0('Deenick/',gene,'.perc.expression.boxplot.pdf'))
    print(ggplot(perc.expression, aes(x=celltype, y=NumCellsExpressing, fill=condition)) +
        geom_boxplot() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),
        plot.margin = unit(c(1,1,1,3), "lines")) +
        labs(y='Percentage of cells', x='') + ggtitle(gene))
    dev.off()
}

# Plot distribution of gene for each cell type split by condition
for(gene in surface.features){
    gene_expression <- FetchData(pbmc, vars = gene)
    gene_expression$individual <- pbmc$individual
    gene_expression$condition <- pbmc$condition
    gene_expression$celltype <- pbmc$cellTypist
    gene_sym <- sym(gene)
    pdf(paste0('Deenick/',gene,'.density.pdf'))
    print(ggplot(gene_expression, aes(x=!!gene_sym, fill=condition)) +
        geom_density(alpha=0.5) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),
        plot.margin = unit(c(1,1,1,3), "lines")) +
        labs(y='Density', x='') + ggtitle(gene) +
        facet_wrap(~celltype))
    dev.off()
}


# Correlation of features 
for(cell in names(features)){
    pbmc.markers <- subset(pbmc, idents=cell, features=features[[cell]]$Features)
    psuedobulk <- AggregateExpression(pbmc.markers, group.by='individual')$RNA
    psuedobulk <- scale(psuedobulk)
    colnames(psuedobulk) <- gsub('data\\[, 1\\]', '', names(colnames(psuedobulk)))

    # Reorder columns by meta$individual
    psuedobulk <- psuedobulk[match(meta$individual, rownames(psuedobulk)),]

    pdf(paste0('Deenick/', gsub('/| ', '.', cell), '.heatmap.pdf'))
    column_ha = HeatmapAnnotation(condition=meta$condition, col=list(condition=c('control'='blue', 'disease'='red')))
    print(Heatmap(psuedobulk, name='z-score', column_title=cell, show_row_names=F, show_column_names=F))
    dev.off()

    # Correlation
    cor.results <- cor(psuedobulk, method='spearman')

    # Reorder correlation matrix
    cor.results <- cor.results[order(meta$individual),order(meta$individual)]

    cor.results <- cor.results[match(rownames(cor.results), meta$individual),match(colnames(cor.results), meta$individual)]
    pdf(paste0('Deenick/', gsub('/| ', '.', cell), '.heatmap.pdf'))
    column_ha = HeatmapAnnotation(condition=meta$condition, col=list(condition=c('control'='blue', 'disease'='red')))
    print(Heatmap(psuedobulk, name='z-score', column_title=cell, show_row_names=F, show_column_names=F,
    col=colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red')),
    clustering_distance_columns='euclidean', clustering_method_columns='ward.D2',
    clustering_distance_rows='euclidean', clustering_method_rows='ward.D2'))
    dev.off()
}

pbmc <- ScaleData(pbmc, features=features[[cell]]$Features)
pdf('Deenick/test.heatmap.pdf')
DoHeatmap(pbmc, group.by='condition', features=features[[cell]]$Features)
dev.off()

# PCA of samples
for(cell in names(features)){
    pbmc.markers <- subset(pbmc, idents=cell, features=features[[cell]]$Features)
    psuedobulk <- AggregateExpression(pbmc.markers, group.by='individual')$RNA
    psuedobulk <- scale(psuedobulk)
    meta <- unique(pbmc.markers@meta.data[,c('individual', 'condition')])
    # PCA
    pca <- prcomp(t(psuedobulk), center=F, scale=F)
    pdf(paste0('Deenick/', gsub('/| ', '.', cell), '.pca.pdf'))
    print(ggplot(data.frame(pca$x), aes(x=PC1, y=PC2)) +
        geom_point() +
        labs(x=paste0('PC1 (', round(pca$sdev[1], 2), '%)'), y=paste0('PC2 (', round(pca$sdev[2], 2), '%)')) +
        theme_bw())
    dev.off()  
}

pbmc.subset <- subset(pbmc, condition=='disease', idents='Classical monocytes', features=features[['Classical monocytes']]$Features)
pbmc.subset <- ScaleData(pbmc.subset)
pca.subset <- RunPCA(pbmc.subset, features=features[['Classical monocytes']]$Features)

# Visualize the results
pdf('PCA_monocytes.pdf')
DimPlot(pca_monocytes, group.by = "ident", label = TRUE, repel = TRUE)
dev.off()


library(enrichR)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("Enrichr") # Human genes   
}
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2023", "MSigDB_Hallmark_2020", "KEGG_2021_Human")

enrichR.results <- lapply(features, function(x){
    enriched <- enrichr(x$Features, dbs)
    enriched <- lapply(enriched, function(y){
        subset(y, Adjusted.P.value < 0.05)
    })
})

lapply(enrichR.results, function(x){
    lapply(x, function(y){
        if(nrow(y) > 0){
            plotEnrich(y, showTerms = nrow(y), numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
        }
    })
})

pdf('Naive.B.cells.GPBP.pdf')
plotEnrich(enrichR.results$perm.Naive.B.cells.chrX.txt$GO_Biological_Process_2023, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
dev.off()

pdf('Tcm.Naive.cytotoxic.T.cells.GPBP.pdf')
plotEnrich(enrichR.results$perm.Tcm.Naive.cytotoxic.T.cells.chrX.txt$GO_Biological_Process_2023, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
dev.off()

library(ggplot2)
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')

chrX.files <- list.files('psuedobulk/features/', pattern='enet', full.names=TRUE)
chrX.files <- grep('chrX', chrX.files, value=TRUE)

chrX.coef <- lapply(chrX.files, read.csv, row.names=1)
names(chrX.coef) <- replace.names(gsub('psuedobulk/features//enet_features\\.|\\.chrX.csv', '', chrX.files))
chrX.coef <- chrX.coef[!is.na(names(chrX.coef))]

# bind rows
combined.coef <- dplyr::bind_rows(chrX.coef, .id='celltype')

library(ggrepel)

pdf('APR/feature.coef.violin.pdf')
ggplot(combined.coef, aes(x=celltype, y=coef)) +
    geom_violin() +
    geom_text_repel(data = subset(combined.coef, abs(coef) > 0.05), aes(label = gsub('\\.+\\d+', '', rownames(subset(combined.coef, abs(coef) > 0.05))), x = celltype, y = coef), size = 3, hjust = 1, vjust = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),
    plot.margin = unit(c(1,1,1,2), "lines")) +
    labs(y='Elastic Net Coefficient', x='')
dev.off()

pdf('APR/feature.coef.density.pdf')
ggplot(combined.coef, aes(x=coef)) +
  geom_density() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),
  plot.margin = unit(c(1,1,1,2), "lines")) +
  labs(y='Density', x='Elastic Net Coefficient')
dev.off()

# Adding DEG results to selected features
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
deg <- deg.list('differential.expression/edgeR', filter=F)
names(deg) <- replace.names(gsub('_', '.', names(deg)))

perm.features.file <- list.files('psuedobulk/ML.models/ensemble/features/', pattern='perm', full.names=TRUE)
perm.features.file <- grep('chrX', perm.features.file, value=TRUE)
# select for high performing models
perm.features.file <- perm.features.file[c(5,11,18,21,22)]
perm.features <- lapply(perm.features.file, read.delim, header=T)
names(perm.features) <- replace.names(gsub('perm.|.chrX.txt', '', basename(perm.features.file)))

# Subset DEG results to selected features
deg.features <- lapply(names(perm.features), function(x){
    subset(deg[[x]], gene %in% perm.features[[x]]$Features)
})
names(deg.features) <- names(perm.features)

# Write DEG results to file
lapply(names(deg.features), function(x){
    write.table(deg.features[[x]], paste0('Deenick/', gsub('-|/| ', '.', x), '.chrX.txt'), sep='\t', row.names=F, quote=F)
})

# Plot volcano plots of DEG results
library(ggrepel)
lapply(names(deg.features), function(x){
    pdf(paste0('Deenick/',gsub('-|/| ', '.', x), '.volcano.chrX.pdf'))
    print(ggplot(deg.features[[x]], aes(x=logFC.disease_vs_control, y=-log10(FDR))) +
        geom_point(aes(color=ifelse(FDR < 0.05 & abs(logFC.disease_vs_control) > 0.2, 'Significant', 'Non-significant'))) +
        geom_text_repel(data=subset(deg.features[[x]], FDR < 0.05 & abs(logFC.disease_vs_control) > 0.2), 
        aes(label=gene), size=3, vjust=1, box.padding = 0.5, point.padding = 0.5) +
        labs(color='FDR < 0.05 & | logFC | > 0.2') +
        labs(x='log2 Fold Change', y='-log10 FDR') + ggtitle(x) +
        geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "black", alpha=1) +
        geom_vline(xintercept = -0.2, linetype="dashed", color = "black", alpha=1) +
        geom_vline(xintercept = 0.2, linetype="dashed", color = "black", alpha=1))
    dev.off()
})

library(Seurat)
pbmc <- readRDS('pbmc.female.control-managed.RDS')

pbmc.cell <- subset(pbmc, cellTypist=='Classical monocytes')
# Keep genes with expression in 5% of cells
keep <- rowSums(pbmc.cell@assays$RNA@counts > 0) > ncol(pbmc.cell) * 0.05
features <- names(keep[keep == T])
pbmc.cell <- subset(pbmc.cell, features=features)

cM <- read.delim('psuedobulk/ML.models/ensemble/features/perm.Classical.monocytes.chrX.txt', header=T)

#### Plotting expression of selected features
