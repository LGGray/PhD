# Function to calculate enrichment of selected features for XCI escape genes
library(dplyr)
library(ggplot2)
# Read in ML metrics file
metrics <- read.delim('exp.matrix/metrics/Metrics.combined.txt')

# Filter for F1 > 0.8
metrics.flt <- metrics %>% 
    # filter(F1 >= 0.8) %>%
    mutate(celltype = gsub('.+_', '', model)) %>%
    mutate(features = gsub('^.*\\.', '', model)) %>%
    arrange(celltype)

subset(metrics.flt, features == 'HVG')

# Calculate mean F1 score across celltype
avg.F1 <- metrics.flt %>%
    group_by(celltype) %>%
    summarise(meanF1=mean(F1)) %>%
    mutate(features = gsub('^.*\\.', '', celltype)) %>%
    mutate(celltype = gsub('.HVG|.chrX', '', celltype)) %>%
    data.frame()

# Plot mean F1 score for each feature set
pdf('ML.plots/F1.mean.barplot.all.pdf')
ggplot(avg.F1, aes(x=celltype, y=meanF1, fill=features)) + 
    geom_bar(position="dodge", stat="identity") +
    #rotate plot horizontally
    coord_flip() +
    labs(x='Features', y='Mean F1 score', title='Mean F1 score for each feature set')
dev.off()


# Read in feature files
feature.files <- list.files('ML.models/features/', pattern='.txt', full.names=TRUE)
feature.list <- lapply(feature.files, read.delim)
names(feature.list) <- gsub('_model|.txt', '', basename(feature.files))
feature.list <- feature.list[names(feature.list) %in% metrics.flt$model]



# Calculate enrichment of XCI escape genes
load('../../datasets/XCI/chrX.Rdata')
load('../../datasets/XCI/escapees.Rdata')

enrichment <- lapply(feature.list, function(x){
    a <- sum(x$Features %in% rownames(escape))
    b <- sum(!(x$Features %in% rownames(escape)))
    c <- sum(rownames(chrX)[!(rownames(chrX) %in% x$Features)] %in% rownames(escape))
    d <- sum(!(rownames(chrX)[!(rownames(chrX) %in% x$Features)] %in% rownames(escape)))
    tmp <- chisq.test(matrix(c(a,b,c,d), nrow=2)) 
})
enrichment

# Create table of feature.list size and enrichment p.value
tmp <- lapply(1:length(enrichment), function(x){
    tmp <- data.frame(nFeatures=nrow(feature.list[[x]]), p.value=enrichment[[x]]$p.value, F1=metrics.flt$F1[x], celltype=metrics.flt$celltype[x])
})
names(tmp) <- names(enrichment)
result <- bind_rows(tmp, .id='model')
result$model <- gsub('_.+', '', result$model)

write.table(result, 'exp.matrix/metrics/Metrics.enrichment.txt', row.names=FALSE, quote=F, sep='\t')

# Plot scatter plot
pdf('Xcape.enrichment.scatterplot.pdf')
ggplot(result, aes(x=F1, y=-log10(p.value), size=nFeatures, colour=celltype, shape=as.factor(model))) + 
  geom_point() +
  labs(x='F1', y='-log10(p.value)', size='nFeatures')
dev.off()

 [1] "Age.associated.B.cells.HVG"       "Age.associated.B.cells.chrX"     
 [3] "B.cells.HVG"                      "B.cells.chrX"                    
 [5] "Classical.monocytes.chrX"         "Cycling.T.cells.HVG"             
 [7] "Cycling.T.cells.chrX"             "Non.classical.monocytes.chrX"    
 [9] "Plasma.cells.chrX"                "Plasmablasts.chrX"               
[11] "Regulatory.T.cells.chrX"          "Tem.Temra.cytotoxic.T.cells.chrX"
[13] "Tem.Trm.cytotoxic.T.cells.chrX"