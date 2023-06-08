# Function to calculate enrichment of selected features for XCI escape genes
library(dplyr)
library(ggplot2)
library(reshape2)

# Create ML.plots directory
dir.create('ML.plots')

# Read in ML metrics file
metrics <- read.delim('exp.matrix/metrics/Metrics.combined.txt')

# Filter for F1 > 0.8
metrics.flt <- metrics %>% 
    filter(F1 >= 0.8) %>%
    mutate(celltype = gsub('.+_|.chrX|.HVG|.HVG-.+', '', model)) %>%
    mutate(features = gsub('^.*\\.', '', model)) %>%
    arrange(celltype)

# # Calculate mean F1 score across celltype
# avg.F1 <- metrics.flt %>%
#     filter(!features %in% c('HVG-X', 'HVG-random')) %>%
#     group_by(celltype) %>%
#     summarise(meanF1=mean(F1)) %>%
#     mutate(features = gsub('^.*\\.', '', celltype)) %>%
#     mutate(celltype = gsub('.HVG|.chrX', '', celltype)) %>%
#     data.frame()
# # Plot mean F1 score for each feature set
# pdf('ML.plots/F1.mean.barplot.all.pdf')
# ggplot(avg.F1, aes(x=celltype, y=meanF1, fill=features)) + 
#     geom_bar(position="dodge", stat="identity") +
#     #rotate plot horizontally
#     coord_flip() +
#     labs(x='Features', y='Mean F1 score', title='Mean F1 score for each feature set')
# dev.off()

# Plot F1 score for each feature set
tmp <- metrics.flt %>%
    filter(!features %in% c('HVG-X', 'HVG-random')) %>%
    mutate(celltype = gsub('.+_|.chrX|.HVG', '', model))
pdf('ML.plots/F1.boxplot.all.pdf')
ggplot(tmp, aes(x=celltype, y=F1, fill=features)) + 
    geom_boxplot() +
    #rotate plot horizontally
    coord_flip() +
    labs(x='Features', y='F1 score', title='F1 scores for trained models')
dev.off()


# Calculate mean F1 score across HVG sets
cells <- unique(metrics.flt[metrics.flt$features %in% c('HVG-X', 'HVG-random'), 'celltype'])
HVG.F1 <- metrics.flt %>%
    filter(grepl('HVG', features)) %>%
    filter(celltype %in% cells) %>%
    data.frame()
# Plot data
pdf('ML.plots/HVG.boxplot.pdf')
ggplot(HVG.F1, aes(x=features, y=F1, fill=features)) + 
    geom_boxplot() +
    geom_jitter() +
    labs(x='Features', y='F1 score', title='F1 scores for HVG sets')
dev.off()

# Plot distributino of F1 scores for each feature set
pdf('ML.plots/HVG.density.pdf')
ggplot(HVG.F1, aes(x=F1, fill=factor(features))) + 
    geom_density(alpha=0.3) +
    labs(x='F1 score', y='Density', title='Density plots for HVG sets')
dev.off()

# Statistical test between HVG and HVG-X, HVG and HVG-random
HVG.F1 %>%
    filter(features %in% c('HVG', 'HVG-X')) %>%
    wilcox.test(F1 ~ features, paired=F, data=.)
HVG.F1 %>%
    filter(features %in% c('HVG', 'HVG-random')) %>%
    wilcox.test(F1 ~ features, paired=F, data=.)


# Read in feature files
feature.files <- list.files('ML.models/features/', pattern='.txt', full.names=TRUE)
feature.list <- lapply(feature.files, read.delim)
names(feature.list) <- gsub('_model|.txt', '', basename(feature.files))
feature.list <- feature.list[names(feature.list) %in% metrics.flt$model]

# Calculate enrichment of XCI escape genes
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')

feature.list <- feature.list[grep('HVG-X|HVG-random', names(feature.list), invert=TRUE)]

enrichment <- lapply(names(feature.list), function(x){
    tmp <- readRDS(paste0('exp.matrix/', gsub('.+_', '', x), '.RDS'))
    features <- feature.list[[x]]$Features
    background <- colnames(tmp)[-c(1:2)]
    background <- background[!(background %in% features)]
    a <- sum(features %in% rownames(escape))
    c <- sum(!(features %in% rownames(escape)))
    b <- sum(background %in% rownames(escape))
    d <- sum(!(background %in% rownames(escape)))
    chisq.test(matrix(c(a,b,c,d), nrow=2)) 
})
names(enrichment) <- names(feature.list)

# enrichment <- enrichment[!sapply(enrichment, function(x) is.na(x$p.value))]
enrich_sig <- enrichment[sapply(enrichment, function(x) x$p.value < 0.05)]
# Create table of feature.list size and enrichment p.value
tmp <- lapply(names(enrich_sig), function(x){
    tmp <- data.frame(nFeatures=nrow(feature.list[[x]]), p.value=enrich_sig[[x]]$p.value)
})
names(tmp) <- names(enrich_sig)
result <- bind_rows(tmp, .id='file') 
result$feature <- gsub('.+\\.', '', result$file)
result$celltype <- gsub('^.+_|.chrX|.HVG', '', result$file)

write.table(result, 'exp.matrix/metrics/Metrics.enrichment.txt', row.names=FALSE, quote=F, sep='\t')

avg.F1 <- metrics.flt %>%
    filter(!features %in% c('HVG-X', 'HVG-random')) %>%
    group_by(celltype, features) %>%
    summarise(meanF1=mean(F1)) %>%
    mutate(file=paste0(celltype, '.', features)) %>%
    data.frame()
avg.F1 <- merge(result, avg.F1, by='file', all.x=TRUE)
# Plot scatter plot
pdf('ML.plots/Xcape.enrichment.scatterplot.pdf')
ggplot(avg.F1, aes(x=meanF1, y=-log10(p.value), size=nFeatures, colour=factor(celltype.x), shape=factor(features))) + 
  geom_point() +
  labs(x='F1', y='-log10(p.value)', size='nFeatures')
dev.off()

# Create heatmap of selected features
# Read in feature files
feature.files <- list.files('ML.models/features/', pattern='.txt', full.names=TRUE)
feature.list <- lapply(feature.files, read.delim)
names(feature.list) <- gsub('_model|.txt', '', basename(feature.files))
feature.list <- feature.list[names(feature.list) %in% metrics.flt$model]

# Remove duplicated lists
names(feature.list) <- gsub('.+_', '', names(feature.list))
feature.list <- feature.list[!duplicated(names(feature.list))]
features <- unique(unlist(feature.list))
# Select chrX features
features <- features[features %in% rownames(chrX)]
# Create matrix
mtx <- matrix(0, nrow=length(features), ncol=length(feature.list))
rownames(mtx) <- features
colnames(mtx) <- names(feature.list)
# Fill matrix
for(i in 1:length(feature.list)){
    mtx[rownames(mtx) %in% feature.list[[i]]$Features, i] <- 1
}
# Plot heatmap
pdf('ML.plots/feature.heatmap.pdf', width=10, height=10)
ggplot(melt(mtx), aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile() +
    scale_fill_gradientn(colors = c("white", "red"), values = c(0, 1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x='Features', y='Models', title='Selected features')
dev.off()

pdf('ML.plots/feature.heatmap.2.pdf')
gplots::heatmap.2(mtx, col = c("white", "red"), notecol = "black", trace = "none", margins = c(11, 5), cexCol = 1, srtCol = 45)
dev.off()

# Gsea of selected features
XIST.dep <- read.delim('../../datasets/XCI/XIST.dependent.txt')
XISR.ind <- read.delim('../../datasets/XCI/XIST.independent.txt')
features[features %in% XIST.dep$XIST.dependent]
features[features %in% XISR.ind$XIST.independent]


# Plot PCA of celltype and features
# Read in Features
features <- read.delim('ML.models/features/RF_model_Naive.B.cells.chrX.txt')$Features
# Subset for selected features and B cells
pbmc.subset <- subset(pbmc, cellTypist=='Naive B cells', features = features)

Idents(pbmc.subset) <- 'condition'
pbmc.subset <- ScaleData(pbmc.subset, features = features)
pbmc.subset <- RunPCA(pbmc.subset, features = features, npcs = 2)
pdf('ML.plots/Naive.B.cells.pca.pdf', width=10, height=10)
DimPlot(pbmc.subset, reduction='pca', label=FALSE)
dev.off()
