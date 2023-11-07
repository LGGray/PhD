library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

### enrichR ###
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')
feature.files <- list.files('psuedobulk/ML.models/ensemble/features/', pattern='perm.*.chrX.txt', full.names=TRUE)
feature.files <- feature.files[c(2,5,11,18,21,22)]
features <- lapply(feature.files, read.delim, header=T)
names(features) <- replace.names(gsub('perm.|.chrX.txt', '', basename(feature.files)))

pbmc <- readRDS('pbmc.female.control-managed.RDS')

surface.features <- c('MSN', 'LAMP2', 'ATP6AP2', 'CXCR3', 'CD40LG', 'IL2RG')

for(cell in names(features)[-1]){
    pbmc.markers <- subset(pbmc, idents=cell, features=features[[cell]]$Features)
    psuedobulk <- AggregateExpression(pbmc.markers, group.by='individual')$RNA
    meta <- unique(pbmc.markers@meta.data[,c('individual', 'condition')])
    # Correlation
    cor.results <- cor(psuedobulk, method='spearman')
    pdf(paste0('Deenick/', gsub('/| ', '.', cell), '.correlation.heatmap.pdf'))
    column_ha = HeatmapAnnotation(condition=meta$condition, col=list(control='blue', disease='red'))
    print(Heatmap(cor.results, name='Spearman correlation',column_title=cell, show_row_names=F, show_column_names=F,
    col=colorRamp2(c(-1,0,1), c('blue','white','red')), top_annotation=column_ha))
    dev.off()   
    #PCA
    pca <- prcomp(t(psuedobulk), scale=T)
    pdf(paste0('Deenick/', gsub('/| ', '.', cell), '.PCA.pdf'))
    ggplot(data.frame(pca$x), aes(x=PC1, y=PC2, color=factor(meta$condition))) + 
    geom_point() + theme_bw() + ggtitle(cell) + labs(color = "Condition")
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

