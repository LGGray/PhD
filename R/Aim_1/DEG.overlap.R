library(ComplexHeatmap)
library(circlize)

source('../PhD/functions/edgeR.list.R')
load('../datasets/XCI/chrX.Rdata')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')

pSS <- deg.list('pSS_GSE157278/differential.expression/edgeR/', logfc=0.5)
UC <- deg.list('UC_GSE125527/differential.expression/edgeR/', logfc=0.5)
CD_colon <- deg.list('CD_Kong/colon/differential.expression/edgeR/', logfc=0.5)
CD_TI <- deg.list('CD_Kong/TI/differential.expression/edgeR/', logfc=0.5)
SLE <- deg.list('lupus_Chun/differential.expression/edgeR/', logfc=0.5)

# Plot the number of rows in each list
pdf('DEG.overlap.pdf')
barplot(unlist(lapply(list(pSS, UC, CD_colon, CD_TI, SLE.2, SLE.5), function(x) length(x))), 
        names.arg=c('pSS', 'UC', 'CD colon', 'CD TI', 'SLE 0.2', 'SLE 0.5'), 
        xlab='Disease', ylab='Number of DEGs', las=2, cex.names=0.8)
dev.off()

# Barplot of the number of DEGS in SLE.2 and SLE.5 coloured by celltype



plot.data <- data.frame(sle.2=unlist(lapply(SLE.2, nrow)), sle.5=unlist(lapply(SLE.5, nrow))) %>%
  reshape2::melt() %>%
  mutate(variable = factor(variable, levels=c('sle.2', 'sle.5')))

pdf('SLE.2.vs.SLE.5.pdf')
ggplot(plot.data, aes(x=variable, y=value)) + 
geom_bar(stat='identity') + 
theme_bw() + 
theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
ylab('Number of DEGs') + xlab('Cell type') + ggtitle('SLE: logFC > 0.2 vs logFC > 0.5') +
theme(plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# Find common celltypes as names of the lists
common <- Reduce(intersect, list(names(pSS), names(UC), names(CD_colon), names(CD_TI), names(SLE)))
common <- replace.names(gsub('_', '.', gsub('CD16__NK_cells', 'CD16-_NK_cells', common)))

pSS <- pSS[common]
UC <- UC[common]
CD_colon <- CD_colon[common]
CD_TI <- CD_TI[common]
SLE <- SLE[common]

studies <- list('pSS'=pSS, 'UC'=UC, 'CD_colon'=CD_colon, 'CD_TI'=CD_TI, 'SLE'=SLE)

for(study in names(studies)){
  for(celltype in names(studies[[study]])){
    print('C4' %in% studies[study][[celltype]]$gene)
  }
}

disco.score <- function(fc1, fc2, p1, p2) {
  # Calculate the directionality factor, which is -1 if the genes are discordant
  directionality <- ifelse(sign(fc1) != sign(fc2), -1, 1)
  # Calculate the disco.score with the directionality factor included
  return(directionality * log2(abs(fc1)) * log2(abs(fc2)) * abs(log10(p1) + log10(p2)))
}

# Combinations of studies
study_combinations <- combn(names(studies), 2, simplify = FALSE)

disco.results <- list()
for( cell in common){
  # Initialize an empty list to store results
  results_list <- list()

  # Iterate over each combination of studies
  for (i in seq_along(study_combinations)) {
    # Extract the study names for this combination
    study_names <- study_combinations[[i]]
    
    # Merge the corresponding data frames by 'gene'
    merged <- merge(studies[[study_names[1]]][[cell]], studies[[study_names[2]]][[cell]], by = 'gene', all = FALSE)
    
    # Calculate the disco.score for each gene
    result <- data.frame(gene = merged$gene, 
                        disco.score = disco.score(merged$logFC.disease_vs_control.x, 
                                                  merged$logFC.disease_vs_control.y, 
                                                  merged$FDR.x, 
                                                  merged$FDR.y))
    
    # Store the result in the list with a name that reflects the combination
    results_list[[paste(study_names, collapse = "_")]] <- result
  }
}
plot.data <- dplyr::bind_rows(results_list, .id='combination')
plot.data_wide <- reshape2::dcast(plot.data, gene ~ combination)
rownames(plot.data_wide) <- plot.data_wide$gene
plot.data_wide <- plot.data_wide[,-1]
plot.data_wide[is.na(plot.data_wide)] <- 0

pdf('test.disco.heatmap.pdf')
Heatmap(as.matrix(plot.data_wide), show_row_names=FALSE)
dev.off()


exp <- readRDS('psuedobulk/Naive.B.cells.chrX.RDS')
class <- exp$class
features <- read.delim('psuedobulk/ML.models/ensemble/features/perm.Naive.B.cells.chrX.txt')
exp <- exp[, features[features$Features %in% rownames(escape),]]
expression.matrix <- matrix(rpois(1000, 10), nrow=100, ncol=10)

pdf('test.heatmap.pdf')
column_ha = HeatmapAnnotation(bar=class)
Heatmap(t(scale(exp)), clustering_distance_columns = "spearman", clustering_distance_rows = "spearman", 
        clustering_method_columns = "complete", clustering_method_rows = "complete", 
        show_row_names = FALSE, show_column_names = FALSE, name='logFC z-score',
        top_annotation = column_ha)
dev.off()



ggplot(plot.data, aes(x=gene, y=combination, fill=disco.score)) + geom_tile()
dev.off()


# Create heatmap of chrX expression in each celltype and disease
common_genes <- lapply(common, function(celltype) {
  unique(unlist(list(pSS[[celltype]]$gene, UC[[celltype]]$gene, SLE[[celltype]]$gene, CD_colon[[celltype]]$gene, CD_TI[[celltype]]$gene)))
})
names(common_genes) <- common

# Upset plot for gene overlap
library(UpSetR)
upset.list <- fromList(common_genes)
rownames(upset.list) <- unique(unlist(common_genes))
pdf('chrX.upregulated.upset.pdf')
upset(upset.list, nsets=8, nintersects=NA)
dev.off()

upset.list$count <- rowSums(upset.list)
rownames(upset.list)[order(upset.list$count, decreasing=TRUE)][1:10]





matrices <- lapply(common, function(celltype) {
    genes <- common_genes[[celltype]]
    mtx <- matrix(0, nrow=length(common_genes[[celltype]]), ncol=5)
    rownames(mtx) <- genes
    colnames(mtx) <- c('pSS', 'UC', 'CD colon', 'CD TI', 'SLE')
    mtx[,1] <- pSS[[celltype]]$logFC.disease_vs_control[match(genes, pSS[[celltype]]$gene)]
    mtx[,2] <- UC[[celltype]]$logFC.disease_vs_control[match(genes, UC[[celltype]]$gene)]
    mtx[,3] <- SLE[[celltype]]$logFC.disease_vs_control[match(genes, SLE[[celltype]]$gene)]
    mtx[,4] <- CD_colon[[celltype]]$logFC.disease_vs_control[match(genes, CD_colon[[celltype]]$gene)]
    mtx[,5] <- CD_TI[[celltype]]$logFC.disease_vs_control[match(genes, CD_TI[[celltype]]$gene)]
    mtx[is.na(mtx)] <- 0
    return(mtx)
})
names(matrices) <- common

# Correlation of logFC
library(reshape2)
library(ggplot2)

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

reorder_cormat <- function(cormat){
# Use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
}

[1] "DC2"                         "Memory_B_cells"             
[3] "Naive_B_cells"               "Plasma_cells"               
[5] "Regulatory_T_cells"          "Tcm_Naive_helper_T_cells"   
[7] "Tem_Effector_helper_T_cells" "Tem_Trm_cytotoxic_T_cells" 

cormat <- cor(matrices[[8]], method='spearman')
cormat[is.na(cormat)] <- 0
# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
pdf('Tem_Trm_cytotoxic_T_cells.correlation.pdf')
ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
    name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
    ggtitle("Tem/Trm cytotoxic T cells")+ xlab("") + ylab("") +
 coord_fixed()
dev.off()

pdf('DC2.correlation.pdf')
Heatmap(correlation, upper='triangle')
dev.off()

pdf('Memory_B_cells.chrX.heatmap.pdf', width=10, height=11)
Heatmap(scale(matrices[[2]]), clustering_distance_rows = "spearman", clustering_distance_columns = "spearman",
clustering_method_rows = "average", clustering_method_columns = "average",
col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), name='logFC z-score', 
column_title = "Memory B cells", column_title_side = "top",
column_names_rot = 45, column_names_side = "bottom", column_dend_side = "bottom", show_row_names = TRUE)
dev.off()


pdf('CD2.chrX.heatmap.pdf', width=10, height=11)
ggplot(melt(matrices[[7]]), aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low='blue', high='red') + 
  theme_bw() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  ylab("Gene") + xlab("Cell type") + ggtitle("Tem/Trm cytotoxic T cells: chrX logFC ") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

library(gplots)

# Define a function to create a heatmap for a given matrix
create_heatmap <- function(mtx) {
  heatmap.2(mtx, trace='none', col=rev(colorRampPalette(c('blue', 'white', 'red'))(100)), Rowv=FALSE, Colv=FALSE, margins=c(10,10))
}

pdf('DC2.chrX.heatmap.pdf')
heatmap.2(matrices[[1]], trace='none', col=rev(colorRampPalette(c("blue", "white", "red"))(100)), 
          scale='none', margins=c(10,10), main='DC2 chrX genes', xlab='', ylab='',
          srtCol=45, key=FALSE, dendrogram='none')
dev.off()

# Create a heatmap for each cell type
for (celltype in common) {
  mtx <- matrices[[celltype]]
  create_heatmap(mtx)
}

####################
chrX.mtx <- lapply(1:length(common), function(i){
  pss <- subset(pSS[[common[i]]], gene %in% rownames(chrX))$gene
  uc <- subset(UC[[common[i]]], gene %in% rownames(chrX))$gene
  sle <- subset(SLE[[common[i]]], gene %in% rownames(chrX))$gene
  cd_colon <- subset(CD_colon[[common[i]]], gene %in% rownames(chrX))$gene
  cd_ti <- subset(CD_TI[[common[i]]], gene %in% rownames(chrX))$gene
  lst <- fromList(list(pss, uc, sle, cd_colon, cd_ti))
  rownames(lst) <- unique(unlist(list(pss, uc, sle, cd_colon, cd_ti)))
  return(lst)
})
names(chrX.mtx) <- common

sort(table(unlist(sapply(chrX.mtx, function(x) rownames(x)))))

