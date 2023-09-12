library(gplots)
library(ggplot2)
library(dplyr)
# list files in data.splits but ignore X_test, X_tune, X_train, y_test, y_tune, y_train
files <- list.files('data.splits', pattern='.csv', full.names=TRUE)
files <- files[!grepl('X_test|X_tune|X_train|y_test|y_tune|y_train', files)]
# Remove empty files
files <- files[sapply(files, function(x) file.size(x) > 1)]

chrX.files <- grep('chrX', files, value=TRUE)
disgene.files <- grep('disgene', files, value=TRUE)
HVG.files <- grep('HVG', files, value=TRUE)

chrX.features <- lapply(chrX.files, read.csv)
names(chrX.features) <- gsub('.chrX.csv', '', basename(chrX.files))
HVG.features <- lapply(HVG.files, read.csv)
names(HVG.features) <- gsub('.HVG.csv', '', basename(HVG.files))
disgene.features <- lapply(disgene.files, read.csv)
names(disgene.features) <- gsub('.disgene.csv', '', basename(disgene.files))

# Create matrix of chrX features log2FC
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
degs <- deg.list('differential.expression/edgeR', filter=FALSE)
names(degs) <- gsub('_', '.', names(degs))
degs_chrX.features <- lapply(names(chrX.features), function(x) degs[[x]][match(chrX.features[[x]]$X0, degs[[x]]$gene),'logFC.disease_vs_control'])
chrX.mat <- matrix(0, nrow=length(unique(unlist(chrX.features))), ncol=length(chrX.features))
rownames(chrX.mat) <- unique(unlist(chrX.features))
colnames(chrX.mat) <- names(chrX.features)
for (i in 1:length(chrX.features)) {
  chrX.mat[match(chrX.features[[i]]$X0, rownames(chrX.mat)),i] <- degs_chrX.features[[i]]
}
pdf('APR/chrX.features.heatmap.pdf')
heatmap.2(chrX.mat, trace='none', key=TRUE, col=colorpanel(10, 'blue', 'white', 'red'), 
scale='row', margins=c(10,10), srtCol=45, cexCol=1, cexRow=0.3)
dev.off()

# intersect each list
chrX.disgene <- lapply(names(disgene.features), function(x) intersect(chrX.features[[x]]$X0, disgene.features[[x]]$X0))
names(chrX.disgene) <- names(disgene.features)
HVG.disgene <- lapply(names(disgene.features), function(x) intersect(HVG.features[[x]]$X0, disgene.features[[x]]$X0))
names(HVG.disgene) <- names(disgene.features)
chrX.HVG <- lapply(names(chrX.features), function(x) intersect(chrX.features[[x]]$X0, HVG.features[[x]]$X0))
names(chrX.HVG) <- names(chrX.features)

# Create binary matrix of chrX.disgene log2FC
degs_chrX.disgene <- lapply(names(chrX.disgene), function(x) degs[[x]][match(chrX.disgene[[x]], degs[[x]]$gene),'logFC.disease_vs_control'])
chrX.disgene.mat <- matrix(0, nrow=length(unique(unlist(chrX.disgene))), ncol=length(chrX.disgene))
rownames(chrX.disgene.mat) <- unique(unlist(chrX.disgene))
colnames(chrX.disgene.mat) <- names(chrX.disgene)
for (i in 1:length(chrX.disgene)) {
  chrX.disgene.mat[match(chrX.disgene[[i]], rownames(chrX.disgene.mat)),i] <- degs_chrX.disgene[[i]]
}
pdf('APR/chrX.disgene.features.heatmap.pdf')
heatmap.2(chrX.disgene.mat, trace='none', key=TRUE, col=colorpanel(10, 'blue', 'white', 'red'), 
scale='row', margins=c(12,12), srtCol=45, cexCol=1, cexRow=1)
dev.off()

# Create binary matrix of HVG.disgene log2FC
degs_HVG.disgene <- lapply(names(HVG.disgene), function(x) degs[[x]][match(HVG.disgene[[x]], degs[[x]]$gene),'logFC.disease_vs_control'])
HVG.disgene.mat <- matrix(0, nrow=length(unique(unlist(HVG.disgene))), ncol=length(HVG.disgene))
rownames(HVG.disgene.mat) <- unique(unlist(HVG.disgene))
colnames(HVG.disgene.mat) <- names(HVG.disgene)
for (i in 1:length(HVG.disgene)) {
  HVG.disgene.mat[match(HVG.disgene[[i]], rownames(HVG.disgene.mat)),i] <- degs_HVG.disgene[[i]]
}
pdf('APR/HVG.disgene.features.heatmap.pdf')
heatmap.2(HVG.disgene.mat, trace='none', key=TRUE, col=colorpanel(10, 'blue', 'white', 'red'),
scale='row', margins=c(12,12), srtCol=45, cexCol=1)
dev.off()

# Create binary matrix of chrX.HVG log2FC
degs_chrX.HVG <- lapply(names(chrX.HVG), function(x) degs[[x]][match(chrX.HVG[[x]], degs[[x]]$gene),'logFC.disease_vs_control'])
chrX.HVG.mat <- matrix(0, nrow=length(unique(unlist(chrX.HVG))), ncol=length(chrX.HVG))
rownames(chrX.HVG.mat) <- unique(unlist(chrX.HVG))
colnames(chrX.HVG.mat) <- names(chrX.HVG)
for (i in 1:length(chrX.HVG)) {
  chrX.HVG.mat[match(chrX.HVG[[i]], rownames(chrX.HVG.mat)),i] <- degs_chrX.HVG[[i]]
}
pdf('APR/chrX.HVG.features.heatmap.pdf')
heatmap.2(chrX.HVG.mat, trace='none', key=TRUE, col=colorpanel(10, 'blue', 'white', 'red'),
scale='row', margins=c(12,12), srtCol=45, cexCol=1, cexRow=1)
dev.off()


# Read in metrics files for individual chrX models
library(dplyr)
metrics.files <- list.files('exp.matrix/metrics', pattern=c('_metrics_', 'chrX.csv'), full.names=TRUE)
metrics.chrX <- lapply(metrics.files, read.csv)
names(metrics.chrX) <- gsub('.chrX.csv', '', basename(metrics.files))
metrics.chrX.df <- dplyr::bind_rows(metrics.chrX, .id='model')
head(metrics.chrX.df)

# Pull together metrics for chrX
# Get information about files in directory
file.info <- file.info(list.files('exp.matrix/metrics', full.names = TRUE, pattern = '_metrics_'))
# Convert creation dates to POSIXct format
creation_dates <- as.POSIXct(file.info$ctime, origin = '1970-01-01')
# Select files created after a certain date
metrics.files <- list.files('exp.matrix/metrics', full.names = TRUE, pattern = '_metrics_')[creation_dates > as.POSIXct('2023-09-07')]
metrics.chrX <- lapply(grep('chrX', metrics.files, value=TRUE), read.csv)
names(metrics.chrX) <- gsub('.chrX.csv', '', basename(grep('chrX', metrics.files, value=TRUE)))
metrics.chrX.df <- dplyr::bind_rows(metrics.chrX, .id='model')

# Add columns for celltype, features and ML
metrics.chrX.df <- metrics.chrX.df[,-9] %>% 
    mutate(celltype = gsub('.+_metrics_', '', model)) %>%
    mutate(ML = gsub('_metrics_.+', '', model)) %>%
    arrange(celltype)

# Plot forest plot for each cell type
pdf('APR/F1.forest.chrX.pdf')
ggplot(metrics.chrX.df, aes(x=F1, y=celltype, color = ML)) +
    geom_point(size = 1, position=position_jitter(height = 0.5, seed = 42)) +
    geom_errorbarh(
        aes(xmin = F1_lower, xmax = F1_upper),
        height = 0.2,
        position=position_jitter(height = 0.5, seed = 42)) +
    geom_vline(xintercept = 0.8, linetype = 'dotted') +
    theme_bw() +
    xlab("F1 score") + ylab("Cell type") + ggtitle("F1 score for chrX models") +
    labs(color='ML')
dev.off()

# Pull together metrics for disgene
# Get information about files in directory
file.info <- file.info(list.files('exp.matrix/metrics', full.names = TRUE, pattern = '_metrics_'))
# Convert creation dates to POSIXct format
creation_dates <- as.POSIXct(file.info$ctime, origin = '1970-01-01')
# Select files created after a certain date
metrics.files <- list.files('exp.matrix/metrics', full.names = TRUE, pattern = '_metrics_')[creation_dates > as.POSIXct('2023-09-07')]
metrics.disgene <- lapply(grep('disgene', metrics.files, value=TRUE), read.csv)
names(metrics.disgene) <- gsub('.disgene.csv', '', basename(grep('disgene', metrics.files, value=TRUE)))
metrics.disgene.df <- dplyr::bind_rows(metrics.disgene, .id='model')

# Add columns for celltype, features and ML
metrics.disgene.df <- metrics.disgene.df[,-9] %>% 
    mutate(celltype = gsub('.+_metrics_', '', model)) %>%
    mutate(ML = gsub('_metrics_.+', '', model)) %>%
    arrange(celltype)

# Plot forest plot for each cell type
pdf('APR/F1.forest.disgene.pdf')
ggplot(metrics.disgene.df, aes(x=F1, y=celltype, color = ML)) +
    geom_point(size = 1, position=position_jitter(height = 0.5, seed = 42)) +
    geom_errorbarh(
        aes(xmin = F1_lower, xmax = F1_upper),
        height = 0.2,
        position=position_jitter(height = 0.5, seed = 42)) +
    geom_vline(xintercept = 0.8, linetype = 'dotted') +
    theme_bw() +
    xlab("F1 score") + ylab("Cell type") + ggtitle("F1 score for disgene models") +
    labs(color='ML')
dev.off()

# Pull together metrics for chrX ensemble models
ensemble.files <- list.files('ML.models/ensemble', pattern='metrics.*\\chrX.*\\.csv', full.names=TRUE)
metrics.ensemble <- lapply(ensemble.files, read.csv)
names(metrics.ensemble) <- gsub('metrics_|.chrX.csv', '', basename(ensemble.files))
metrics.ensemble.df <- dplyr::bind_rows(metrics.ensemble, .id='model')

ensembl.features <- list.files('data.splits', pattern='chrX', full.names=TRUE)
ensembl.features <- ensembl.features[!grepl('X_test|X_tune|X_train|y_test|y_tune|y_train', ensembl.features)]
features <- lapply(ensembl.features, read.csv)
names(features) <- gsub('.chrX.csv', '', basename(ensembl.features))
feature_size <- unlist(lapply(features[metrics.ensemble.df$model], nrow))

# Plot forest plot for each cell type
pdf('APR/F1.forest.chrX.ensemble.pdf')
ggplot(metrics.ensemble.df, aes(x=F1, y=model, color=AUC, size=feature_size)) +
    geom_point() + 
    geom_vline(xintercept = 0.8, linetype = 'dotted') +
    scale_color_gradient(low = "blue", high = "red") +
    theme_bw() +
    xlab("F1 score") + ylab("Cell type") + ggtitle("F1 score for ensemble chrX models")
dev.off()

# Pull together metrics for disgene ensemble models
ensemble.files <- list.files('ML.models/ensemble', pattern='metrics.*\\disgene.*\\.csv', full.names=TRUE)
metrics.ensemble <- lapply(ensemble.files, read.csv)
names(metrics.ensemble) <- gsub('metrics_|.disgene.csv', '', basename(ensemble.files))
metrics.ensemble.df <- dplyr::bind_rows(metrics.ensemble, .id='model')

ensembl.features <- list.files('data.splits', pattern='disgene', full.names=TRUE)
ensembl.features <- ensembl.features[!grepl('X_test|X_tune|X_train|y_test|y_tune|y_train', ensembl.features)]
features <- lapply(ensembl.features, read.csv)
names(features) <- gsub('.disgene.csv', '', basename(ensembl.features))
feature_size <- unlist(lapply(features[metrics.ensemble.df$model], nrow))

# Plot forest plot for each cell type
pdf('APR/F1.forest.disgene.ensemble.pdf')
ggplot(metrics.ensemble.df, aes(x=F1, y=model, color=AUC, size=feature_size)) +
    geom_point() + 
    geom_vline(xintercept = 0.8, linetype = 'dotted') +
    scale_color_gradient(low = "blue", high = "red") +
    theme_bw() +
    xlab("F1 score") + ylab("Cell type") + ggtitle("F1 score for ensemble disgene models")
dev.off()

subset(degs$Non.classical.monocytes, gene %in% features$Non.classical.monocytes$X0)
subset(degs$Age.associated.B.cells, gene %in% features$Age.associated.B.cells$X0)

# Read in OneK1K object and subset for Non.classical.monocytes and Age.associated.B.cells
library(Seurat)
onek1k <- readRDS('onek1k.RDS')
Idents(onek1k) <- 'predicted.celltype.l3'
metadata <- readRDS('../datasets/OneK1k/disease.meta.RDS')
rownames(metadata) <- metadata$individual

# Remove NAs from metadata
metadata$X08_Autoimmune_Disease[is.na(metadata$X08_Autoimmune_Disease)] <- 0
metadata$X11_Rheumatoid_arthritis[is.na(metadata$X11_Rheumatoid_arthritis)] <- 0
metadata$X12_UlcerativeColitis[is.na(metadata$X12_UlcerativeColitis)] <- 0
metadata$X14_Autoimmune_Disease_Other[is.na(metadata$X14_Autoimmune_Disease_Other)] <- 0

# Count number of autoimmune diseases and create class label
autoimmune <- apply(metadata[,c('X08_Autoimmune_Disease', 'X11_Rheumatoid_arthritis', 'X12_UlcerativeColitis', 'X14_Autoimmune_Disease_Other')], 1, sum) 
autoimmune <- ifelse(autoimmune > 0, 1, 0)
# Add class label to metadata
onek1k$class <- autoimmune[onek1k$individual]

# Read in model features
features <- read.csv('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/lupus_Chun/data.splits/Non.classical.monocytes.chrX.csv')$X0

# Subset for CD16+ monocytes
ncM <- subset(onek1k, idents='CD16 Mono', features=features, sex=='F')
exp <- GetAssayData(ncM, assay='SCT', slot='counts')
exp <- data.frame(t(as.matrix(exp)))
exp <- cbind(class=ncM$class, exp)

# Save expression matrix
saveRDS(exp, 'OneK1K.exp.matrix/Non.classical.monocytes.chrX.RDS')