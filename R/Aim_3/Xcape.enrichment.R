# Function to calculate enrichment of selected features for XCI escape genes
library(dplyr)
library(ggplot2)
library(reshape2)
library(UpSetR)
library(rstatix)

# Create ML.plots directory
dir.create('ML.plots')

# Read in ML metrics file
metrics <- read.delim('exp.matrix/metrics/Metrics.combined.txt')

# Add columns for celltype, features and ML
metrics <- metrics %>% 
    mutate(celltype = gsub('.+_|.chrX|.HVG|.HVG-.+', '', model)) %>%
    mutate(features = gsub('^.*\\.', '', model)) %>%
    mutate(ML = gsub('_.+', '', model)) %>%
    arrange(celltype)

# Calculate mean F1 score across celltype
avg.F1 <- metrics.flt %>%
    filter(!features %in% c('HVG', 'HVG-X', 'HVG-random')) %>%
    group_by(celltype) %>%
    summarise(meanF1=mean(F1)) %>%
    arrange(meanF1) %>%
    mutate(features = gsub('^.*\\.', '', celltype)) %>%
    mutate(celltype = gsub('.chrX', '', celltype)) %>%
    data.frame()
# # Plot mean F1 score for each feature set
# pdf('ML.plots/F1.mean.barplot.all.pdf')
# ggplot(avg.F1, aes(x=celltype, y=meanF1, fill=features)) + 
#     geom_bar(position="dodge", stat="identity") +
#     #rotate plot horizontally
#     coord_flip() +
#     labs(x='Features', y='Mean F1 score', title='Mean F1 score for each feature set')
# dev.off()

# Plot the F1 scores for each model
plot.data <- subset(metrics, features == 'chrX')
pdf('ML.plots/F1.forest.all.pdf')
ggplot(metrics, aes(x=F1, y=celltype, color = ML)) +
    geom_point(size = 1, position=position_jitter(height = 0.5, seed = 42)) +
    geom_errorbarh(
        aes(xmin = F1_lower, xmax = F1_upper),
        height = 0.2,
        position=position_jitter(height = 0.5, seed = 42)) +
    geom_vline(xintercept = 0.8, linetype = 'dotted') +
    theme_bw() +
    xlab("F1 score") + ylab("Cell type") + ggtitle("F1 score for each cell type and model type") +
    labs(color='Features')
dev.off()

# Plot F1 scores for the high performing models (F1 > 0.8)
metrics.flt <- subset(metrics, F1 > 0.8)
plot.data <- subset(metrics.flt, features == 'HVG')
pdf('ML.plots/F1.forest.HVG.filtered.female.pdf')
ggplot(plot.data, aes(x=F1, y=celltype, color = ML)) +
    geom_point(size = 1, position=position_jitter(height = 0.5, seed = 42)) +
    geom_errorbarh(
        aes(xmin = F1_lower, xmax = F1_upper),
        height = 0.2,
        position=position_jitter(height = 0.5, seed = 42)) +
    theme_bw() +
    xlab("F1 score") + ylab("Cell type") + ggtitle("F1 score for each cell type and model type (F1 > 0.8)") +
    labs(color='Features')
dev.off()

# Plot F1 score for each feature set
pdf('ML.plots/F1.boxplot.all.female.pdf')
ggplot(metrics, aes(x=celltype, y=F1, fill=features)) + 
    geom_boxplot() +
    geom_hline(yintercept = 0.8, linetype = 'dotted') +
    #rotate plot horizontally
    coord_flip() +
    labs(x='', y='F1 score', title='F1 scores for trained models')
dev.off()


# Calculate F1 score across HVG sets
HVG.F1 <- metrics %>%
    filter(grepl('HVG', features)) %>%
    mutate(features = factor(features, levels=c('HVG', 'HVG-X', 'HVG-random'))) %>%
    data.frame()
# Plot data
pdf('ML.plots/HVG.boxplot.pdf')
ggplot(HVG.F1, aes(x=features, y=F1, fill=features)) + 
    geom_boxplot() +
    geom_jitter() +
    geom_hline(yintercept = 0.8, linetype = 'dotted') +
    labs(x='Features', y='F1 score', title='F1 scores for HVG sets')
dev.off()

# Plot distribution of F1 scores for each feature set
pdf('ML.plots/HVG.density.pdf')
ggplot(HVG.F1, aes(x=F1, fill=factor(features))) + 
    geom_density(alpha=0.3) +
    labs(x='F1 score', y='Density', title='Density plots for HVG sets')
dev.off()

# Test for difference in F1 scores between chrX models
chrX.F1 <- subset(metrics, features=='chrX')
lapply(split(chrX.F1, chrX.F1$celltype), function(x){
    x$ML <- factor(x$ML)
    if(length(unique(x$ML)) > 1){
        kruskal.test(F1 ~ ML, data=x)
    }
})
kruskal.test(F1 ~ features, data=subset(metrics))
dunn_test(data=metrics, formula=F1 ~ features, p.adjust.method ='fdr')

pdf('ML.plots/all.boxplot.pdf', height=12)
ggplot(metrics, aes(x=features, y=F1, fill=features)) + 
    geom_boxplot() +
    geom_jitter() +
    geom_hline(yintercept = 0.8, linetype = 'dotted') +
    theme(axis.text.x = element_blank()) +
    labs(x='', y='F1 score', title='F1 scores for different feature sets')
dev.off()

pdf('ML.plots/chrX.boxplot.pdf')
ggplot(chrX.F1, aes(x=ML, y=F1, fill=ML)) + 
    geom_boxplot() +
    geom_jitter() +
    geom_hline(yintercept = 0.8, linetype = 'dotted') +
    labs(x='Features', y='F1 score', title='F1 scores for chrX models')
dev.off()

# test for difference in F1 scores between HVG sets split by cell type
lapply(split(HVG.F1, HVG.F1$celltype), function(x){
    x$features <- factor(x$features, levels=c('HVG', 'HVG-X', 'HVG-random'))
    if(length(unique(x$features)) > 1){
        kruskal.test(F1 ~ features, data=x)
    }
})
kruskal.test(F1 ~ features, data=subset(HVG.F1, celltype=='DC2'))
dunn_test(data=subset(HVG.F1, celltype=='DC2'), formula=F1 ~ features, p.adjust.method ='fdr')
# test for difference in F1 scores between HVG sets across all cell types
kruskal.test(F1 ~ features, data=HVG.F1)

metrics.flt <- subset(metrics, F1 > 0.8)

# Read in feature files
feature.files <- list.files('ML.models/features/', pattern='.txt', full.names=TRUE)
feature.list <- lapply(feature.files, read.delim)
names(feature.list) <- gsub('_model|.txt', '', basename(feature.files))
feature.list <- feature.list[names(feature.list) %in% metrics.flt$model]

# Calculate enrichment of XCI escape genes
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')

feature.list <- feature.list[grep('HVG|HVG-X|HVG-random', names(feature.list), invert=TRUE)]

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
result$F1 <- metrics[match(result$file, metrics$model), 'F1']
result$feature <- gsub('.+\\.', '', result$file)
result$celltype <- gsub('^.+_|.chrX|.HVG', '', result$file)

write.table(result, 'exp.matrix/metrics/Metrics.enrichment.txt', row.names=FALSE, quote=F, sep='\t')

# Plot scatter plot
pdf('ML.plots/Xcape.enrichment.scatterplot.pdf')
ggplot(result, aes(x=F1, y=-log10(p.value), size=nFeatures)) + 
  geom_point() +
  geom_text(aes(label=file), size = 3, vjust="inward",hjust="inward") +
  labs(x='F1', y='-log10(p.value)', size='nFeatures')
dev.off()

# avg.F1 <- metrics.flt %>%
#     filter(!features %in% c('HVG-X', 'HVG-random')) %>%
#     group_by(celltype, features) %>%
#     summarise(meanF1=mean(F1)) %>%
#     mutate(file=paste0(celltype, '.', features)) %>%
#     data.frame()
# avg.F1 <- merge(result, avg.F1, by='file', all.x=TRUE)
# # Plot scatter plot
# pdf('ML.plots/Xcape.enrichment.scatterplot.pdf')
# ggplot(avg.F1, aes(x=meanF1, y=-log10(p.value), size=nFeatures, colour=factor(celltype.x), shape=factor(features))) + 
#   geom_point() +
#   labs(x='F1', y='-log10(p.value)', size='nFeatures')
# dev.off()

# Create heatmap of selected features
# Read in feature files
feature.files <- list.files('ML.models/features/', pattern='chrX.txt', full.names=TRUE)
feature.list <- lapply(feature.files, read.delim)
names(feature.list) <- gsub('_model|.txt', '', basename(feature.files))
feature.list <- feature.list[names(feature.list) %in% metrics.flt$model]
features <- unique(unlist(feature.list))

# names(feature.list) <- gsub('.+_', '', names(feature.list))
# feature.list <- feature.list[!duplicated(names(feature.list))]

# Create matrix
mtx <- matrix(0, nrow=length(features), ncol=length(feature.list))
rownames(mtx) <- features
colnames(mtx) <- names(feature.list)
# Fill matrix
for(i in 1:length(feature.list)){
    mtx[rownames(mtx) %in% feature.list[[i]]$Features, i] <- 1
}
# Plot heatmap
pdf('ML.plots/feature.heatmap.chrX.pdf', width=15, height=10)
ggplot(melt(mtx), aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile() +
    scale_fill_gradientn(colors = c("white", "red"), values = c(0, 1)) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x='Features', y='Models', title='Selected features')
dev.off()

pdf('ML.plots/feature.heatmap.2.pdf', width=10, height=10)
gplots::heatmap.2(mtx, col = c("white", "red"), notecol = "black", trace = "none", margins = c(13, 5), cexRow = 1, srtCol = 45, key = FALSE)
dev.off()

# Gsea of selected features
XIST.dep <- read.delim('../../datasets/XCI/XIST.dependent.txt')
XISR.ind <- read.delim('../../datasets/XCI/XIST.independent.txt')
features[features %in% XIST.dep$XIST.dependent]
features[features %in% XISR.ind$XIST.independent]


# Plot PCA of celltype and features
# Read in Features
features <- read.delim('ML.models/features/RF_model_Age.chrX.txt')$Features
# Subset for selected features and B cells
pbmc <- readRDS('pbmc.female.RDS')
pbmc.subset <- subset(pbmc, cellTypist=='Tem/Temra cytotoxic T cells', features = features)

Idents(pbmc.subset) <- 'condition'
pbmc.subset <- ScaleData(pbmc.subset, features = features)
pbmc.subset <- RunPCA(pbmc.subset, features = features, npcs = 2)
pdf('ML.plots/Naive.B.cells.pca.pdf', width=10, height=10)
DimPlot(pbmc.subset, reduction='pca', label=FALSE)
dev.off()


metrics.data <- data.frame(individual=pbmc$individual, condition=factor(pbmc$condition), 
disease_state=pbmc$disease_state, ancestry=pbmc$ethnicity, row.names=NULL)
metrics.data <- unique(metrics.data)
metrics.data$disease_state <- gsub('^na', 'control', metrics.data$disease_state)

age.data <- data.frame(individual=pbmc$individual, disease_state=factor(pbmc$disease_state),
age=pbmc$development_stage, row.names=NULL)
age.data <- unique(age.data)
age.data$age <- as.numeric(gsub('-.+', '', age.data$age))
age.data$disease_state <- gsub('^na', 'control', age.data$disease_state)


# barplot of metrics.data showing proportion of samples in each condition
data_long <- reshape2::melt(metrics.data, id.vars = c("individual", "condition"))
pdf('ML.plots/metrics.data.barplot.pdf')
ggplot(data_long, aes(x = variable, fill = value)) +
  geom_bar(position = "fill") +
  facet_grid(~ condition)
dev.off()

# barplot of age.data showing proportion of samples in each condition
pdf('ML.plots/age.data.barplot.pdf')
ggplot(age.data, aes(x = disease_state, y = age)) +
  geom_boxplot() +
  geom_jitter() +
  scale_y_continuous(breaks=seq(0,80,10))
dev.off()


# Upset plot of selected chrX features
# Read in feature files
feature.files <- list.files('ML.models/features/', pattern='chrX.txt', full.names=TRUE)
feature.list <- lapply(feature.files, function(x) read.delim(x)$Features)
names(feature.list) <- gsub('_model|.txt', '', basename(feature.files))
feature.list <- feature.list[names(feature.list) %in% metrics.flt$model]

# Calculate enrichment of XCI escape genes
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')

feature.list <- feature.list[grep('HVG|HVG-X|HVG-random', names(feature.list), invert=TRUE)]

# Remove duplicated lists
names(feature.list) <- gsub('.+_', '', names(feature.list))
feature.list <- feature.list[!duplicated(names(feature.list))]

lst <- UpSetR::fromList(feature.list)
rownames(lst) <- unique(unlist(feature.list))
lst <- lst[rownames(lst) %in% rownames(escape),]
# Upset plot with top features
lst <- data.frame(t(lst))
lst <- lst[,order(colSums(lst), decreasing=TRUE)]
pdf('ML.plots/upsetplot.chrX.pdf', width=10, height=10)
upset(data.frame(t(lst)), order.by = "freq", nsets = 49, point.size = 2.5, line.size = 1.5, 
              main.bar.color = "black", sets.bar.color = "black", text.scale = 1.5, 
              matrix.color = "black", shade.color = "black")
dev.off()

lst <- unique(unlist)
features <- lst[lst %in% rownames(escape)]
pdf('ML.plots/escape.heatmap.pdf')
DoHeatmap(pbmc, features=features, group.by='condition')
dev.off()


features <- feature.list[[1]][feature.list[[1]] %in% rownames(escape)]
res$threshold <- ifelse(res$FDR < 0.05, 'red', 'black')
pdf('ML.plots/Cycling.T.cells.volcano.pdf')
# colour points if FDR < 0.05
ggplot(res, aes(x=logFC, y=-log10(FDR))) +
    geom_point(colour=res$threshold) +
    geom_text_repel(data=subset(res, FDR < 0.05), aes(label=gene), size=5, color = 'black') +
    ggtitle("Cycling T cells") +
    labs(x = "logFC", y = "-log10 FDR")
dev.off()
