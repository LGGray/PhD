library(speckle)
library(limma)
library(ggplot2)
library(Seurat)
library(dplyr)

pbmc <- readRDS(commandArgs(trailingOnly = TRUE)[1])

load('/directflow/SCCGGroupShare/projects/lacgra/PhD/R/celltype.colours.RData')

# # Clustered bar plot of cell type proportions
# pdf('celltype_props.pdf')
# g <- plotCellTypeProps(pbmc, clusters=pbmc$cellTypist, sample=pbmc$condition)
# ggplot(g$data, aes(x = Samples, y = Proportions, fill = Clusters)) + 
#     geom_bar(stat = "identity") + scale_fill_manual(values=celltype_colours, name='') +
#     theme(axis.text.x = element_text(hjust = 1)) + ggtitle("pSS Cell Type Proportions") +
#     xlab('') + ylab('')
# dev.off()

# Case control cell proportion volcanoplot
pbmc$condition <- factor(pbmc$condition, levels=c('disease', 'control'))
output.logit <- propeller(clusters=pbmc$cellTypist, sample=pbmc$individual, group=pbmc$condition, transform='logit')
output.asin <- propeller(clusters=pbmc$cellTypist, sample=pbmc$individual, group=pbmc$condition, transform='asin')
write.table(output.logit, 'output.logit.txt')
write.table(output.asin, 'output.asin.txt')

# library(ggrepel)
# pdf('celltype_props_volcano.pdf')
# ggplot(output.logit, aes(x=Tstatistic, y=-log10(FDR), colour=factor(BaselineProp.clusters))) +
#     geom_point() + 
#     scale_colour_manual(values=celltype_colours) +
#     theme(legend.position='none') +
#     xlab('T-statistic') + ylab('-log10(FDR)') +
#     geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "red") +  # Add a dotted line at the significance mark
#     geom_text_repel(aes(label=BaselineProp.clusters, color='black'), vjust=-1) +  # Use ggrepel for labels
#     ggtitle('pSS cell type proportions')
# dev.off()

### Comparing cell type proportions between diseases ###
library(ComplexHeatmap)
library(circlize)
# Read in asin output and plot heatmap of cell type t-statistic
MS <- read.delim('MS_GSE193770/output.logit.txt', sep=' ')
pSS <- read.delim('pSS_GSE157278/output.logit.txt', sep=' ')
UC <- read.delim('UC_GSE125527/output.logit.txt', sep=' ')
SLE <- read.delim('lupus_Chun/output.logit.txt', sep=' ')
CO <- read.delim('CD_Kong/colon/output.logit.txt', sep=' ')
TI <- read.delim('CD_Kong/TI/output.logit.txt', sep=' ')

foo <- lapply(list(MS, pSS, UC, SLE, CO, TI), function(x) {
  subset(x, FDR < 0.05)
})[[4]]

subset(foo, Tstatistic > 0)

pb <- AggregateExpression(pbmc, by='cellTypist')$decontXcounts
sort(pb[grep('^CD16$', rownames(pb)),])

cell.prop.lst <- list(MS=MS[,c(1,6,8)], pSS=pSS[,c(1,6,8)], UC=UC[,c(1,6,8)], SLE=SLE[,c(1,6,8)], CO=CO[,c(1,6,8)], TI=TI[,c(1,6,8)])
# # Replace T statistic with 0 if FDR < 0.05
# cell.prop.lst <- lapply(cell.prop.lst, function(x) {
#   x[x[,3] > 0.05, 2] <- 0
#   x
# })

foo <- dplyr::bind_rows(cell.prop.lst, .id='dataset')
rownames(foo) <- NULL
# if FDR < 0.05 foo$significant == '*'; FDR < 0.01 foo$significant == '**'; FDR < 0.001 foo$significant == '***'
foo$significant <- ifelse(foo$FDR < 0.05, '*', ifelse(foo$FDR < 0.01, '**', ifelse(foo$FDR < 0.001, '***', '')))

celltypes <- unique(unlist(lapply(cell.prop.lst, function(x) rownames(x))))
# Create matrix
cell.prop.mat <- matrix(0, nrow=length(celltypes), ncol=length(cell.prop.lst))
rownames(cell.prop.mat) <- celltypes
colnames(cell.prop.mat) <- names(cell.prop.lst)

# Fill matrix
for (i in 1:length(cell.prop.lst)) {
    cell.prop.mat[match(cell.prop.lst[[i]][,1], rownames(cell.prop.mat)), i] <- cell.prop.lst[[i]][,2]
}

# significance matrix
cell.prop.mat.sig <- matrix(1, nrow=length(celltypes), ncol=length(cell.prop.lst))
rownames(cell.prop.mat.sig) <- celltypes
colnames(cell.prop.mat.sig) <- names(cell.prop.lst)

# Fill matrix
for (i in 1:length(cell.prop.lst)) {
    cell.prop.mat.sig[match(cell.prop.lst[[i]][,1], rownames(cell.prop.mat.sig)), i] <- ifelse(cell.prop.lst[[i]][,3] > 0.05, '', 
    ifelse(cell.prop.lst[[i]][,3] < 0.05 & cell.prop.lst[[i]][,3] > 0.01, '*', 
    ifelse(cell.prop.lst[[i]][,3] < 0.01 & cell.prop.lst[[i]][,3] > 0.001, '**', '***'))) 
}



# Replace FDR with *
cell.prop.mat.sig[cell.prop.mat.sig > 0.05] <- ''
cell.prop.mat.sig[cell.prop.mat.sig < 0.05 & cell.prop.mat.sig > 0.01] <- '*'
cell.prop.mat.sig[cell.prop.mat.sig < 0.01 & cell.prop.mat.sig > 0.001] <- '**'
cell.prop.mat.sig[cell.prop.mat.sig < 0.001] <- '***'

# Plot heatmap
pdf('Aim_1/celltype_props_heatmap.pdf')
col=colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
Heatmap(cell.prop.mat, name='T-statistic', col=col, column_title='', row_title='',
column_names_rot=0, row_names_gp=gpar(fontsize=8))
dev.off()

# Plot heatmap
pdf('Aim_1/celltype_props_heatmap.pdf')
col=colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
Heatmap(cell.prop.mat, name='T-statistic', col=col, column_title='', row_title='',
column_names_rot=0, row_names_gp=gpar(fontsize=8))
dev.off()