library(Seurat)
library(ggplot2)
library(dplyr)

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/chisq.test.degs.R')

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/RA_SDY998/')

# Load the data
pbmc <- readRDS('pbmc.female.RDS')

degs <- edgeR.list('psuedobulk', logfc=0)
names(degs) <- gsub('.edgeR-LRT', '', names(degs))

deg.metrics <- lapply(degs, function(x){
    up.all <- x[which(x$logFC > 0),]
    down.all <- x[which(x$logFC < 0),]
    XCI <- subset(x, gene %in% rownames(chrX))
    up.XCI <- XCI[which(XCI$logFC > 0),]
    down.XCI <- XCI[which(XCI$logFC < 0),]
    data.frame(up=nrow(up.all), down=nrow(down.all), up.chrX=nrow(up.XCI), down.chrX=nrow(down.XCI))
}) %>% do.call(rbind, .) 
deg.metrics

# Get DEGs for chrX
lapply(degs, function(x){
    subset(x, gene %in% rownames(chrX) & logFC > 0)
})
lapply(degs, function(x){
    subset(x, gene %in% rownames(chrX) & logFC < 0)
})

# Count number of cells in each condition for each cell
table(pbmc$condition, pbmc$cellTypist)

ABC <- subset(pbmc, cellTypist == 'Age-associated B cells'))
Bmem <- subset(pbmc, cellTypist == 'Memory B cells')
Treg <- subset(pbmc, cellTypist == 'Regulatory T cells')
Temra <- subset(pbmc, cellTypist == 'Tem/Temra cytotoxic T cells')
Th <- subset(pbmc, cellTypist == 'Tem/Trm cytotoxic T cells' & condition == 'disease')

cell.sum <- rowSums(Temra)
sort(cell.sum[grep('^CD', names(cell.sum))])

# Find cell type markers
pbmc.disease <- subset(pbmc, condition == 'disease')
cell.markers <- FindAllMarkers(pbmc.disease, assay='SCT', min.pct=0.25, logfc.threshold=0.25, only.pos=TRUE)

markers.split <- split(cell.markers, cell.markers$cluster)
lapply(markers.split, function(x){
    x[1:10,7]
})

save(cell.markers, file='disease.cell.markers.RDS')
