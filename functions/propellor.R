library(speckle)
library(limma)
library(ggplot2)
library(Seurat)
library(dplyr)

pbmc <- readRDS(commandArgs(trailingOnly = TRUE)[1])
disease <- commandArgs(trailingOnly = TRUE)[2]

load('/directflow/SCCGGroupShare/projects/lacgra/PhD/R/celltype.colours.RData')
colours <- colours[unique(pbmc$cellTypist)]
                  
pdf('APR/pbmc_celltype_props.pdf')
g <- plotCellTypeProps(pbmc, clusters=pbmc$cellTypist, sample=pbmc$individual, transform='asin')
ggplot(g$data, aes(x = Samples, y = Proportions, fill = Clusters)) + 
    geom_bar(stat = "identity") + scale_fill_manual(values=colours, name='') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("UC Cell Type Proportions") +
    xlab('') + ylab('')
dev.off()

meta <- unique(pbmc@meta.data[,c('individual', 'condition')])
g$data$condition <- meta$condition[match(g$data$Samples, meta$individual)]
g$data$condition <- factor(g$data$condition, levels=c('control', 'disease'))

pdf('APR/celltype_props_boxplot.pdf')
ggplot(g$data, aes(x=Clusters, y=Proportions, fill=Clusters, color=condition)) + 
    geom_boxplot() + scale_fill_manual(values=colours, name='') +
    scale_color_manual(values=c('black', 'black'), name='') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position="none") +
    xlab('') + ylab('')
dev.off()

# Case control cell proportion scatterplot
output.logit <- propeller(clusters=pbmc$cellTypist, sample=pbmc$individual, group=pbmc$condition, transform='logit')
output.asin <- propeller(clusters=pbmc$cellTypist, sample=pbmc$individual, group=pbmc$condition, transform='asin')
pdf('APR/celltype_props_scatterplot.pdf', width=10, height=10)
ggplot(output.logit, aes(x=PropMean.control, y=PropMean.disease, colour=factor(BaselineProp.clusters))) +
    geom_point() + geom_abline(intercept=0, slope=1, linetype="dotted") +
    scale_colour_manual(values=colours, name='') +
    theme(legend.position='bottom', legend.text=element_text(size=8, color="black")) +
    xlab('Control proportion') + ylab('Disease proportion') +
    geom_text(aes(label=ifelse(FDR<0.05, rownames(output.logit), "")), vjust=-1) +
    guides(colour = guide_legend(override.aes = list(shape = 16))) +
    ggtitle('Case-Control Cell Type Proportions')
dev.off()

meta <- data.frame(unique(pbmc@meta.data[,c('individual', 'condition')]))
cell.perc <- pbmc@meta.data %>%
    group_by(individual, cellTypist) %>%
    summarise(count = n()) %>%
    mutate(perc = count / sum(count) * 100) %>%
    left_join(meta, by = "individual")

pdf('APR/celltype_props.violin.asin.pdf', width=10, height=10)
ggplot(cell.perc, aes(x=condition, y=perc, fill=cellTypist)) +
    geom_violin() +
    geom_jitter(width=0.1, size=1, alpha=0.5) +
    geom_boxplot(width=0.1) +
    theme_bw() +
    theme(legend.position = 'none') +
    labs(x='', y='log10 Proportion (%)', size = 14) +
    scale_fill_manual(values=colours[unique(cell.perc$cellTypist)]) +
    scale_y_log10() +
    facet_wrap(~cellTypist)
dev.off()

g + scale_fill_manual(values=colours) +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle(paste("pSS Cell Type Proportions")) +
xlab('') + ylab('') + scale_fill_discrete(name='')
dev.off()

pbmc$condition <- factor(pbmc$condition, levels=c('disease', 'control'))
pbmc$age <- as.numeric(gsub('-year-old human stage', '', pbmc$development_stage))

pbmc@meta.data$condition <- factor(pbmc$condition, levels=c('disease', 'control'))
output.logit <- propeller(clusters=pbmc$cellTypist, sample=pbmc$individual, group=pbmc$condition, transform='logit')
output.asin <- propeller(clusters=pbmc$cellTypist, sample=pbmc$individual, group=pbmc$condition, transform='asin')

props <- getTransformedProps(pbmc$cellTypist, pbmc$individual, transform="asin")
targets = unique(data.frame(condition = pbmc$condition,
                      individual = pbmc$individual,
                      age = pbmc$age))
design <- model.matrix(~ 0 + condition, data=targets)
mycontr <- makeContrasts(conditiondisease-conditioncontrol, levels=design)
result <- propeller.ttest(props, design, contrasts = mycontr, robust=TRUE, trend=FALSE, 
                sort=TRUE)
write.table(result, 'propellor.asin.condition.abundance.txt', sep='\t', quote=FALSE)

# Plot scatterplot of T-statistic vs -log10(FDR). colour and label if significant
pdf('pbmc_celltype_props_diff_abundance.pdf')
ggplot(output.asin, aes(x=Tstatistic, y=-log10(FDR))) + 
    geom_point(aes(colour=ifelse(FDR<0.05, "red", "black"))) +
    geom_text(aes(label=ifelse(FDR<0.05, rownames(output.asin), '')), position = position_jitter(width = 1.5, height = 1.5, ), vjust="inward",hjust="inward") +
    ylab("-log10(FDR)") + xlab("T-statistic") + 
    # rename the legend
    scale_colour_manual(name="Significant", values=c("black", "red"), labels=c("No", "Yes")) +
    ggtitle(paste(disease, " Differential Abundance"))
dev.off()

saveRDS(output.asin, 'propellor.asin.RDS') 


library(ComplexHeatmap)
library(circlize)
# Read in asin output and plot heatmap of cell type t-statistic
pSS <- read.delim('pSS_GSE157278/propeller.asin.txt', sep=' ')
UC <- read.delim('UC_GSE125527/propeller.asin.txt', sep=' ')
SLE <- read.delim('lupus_Chun/propellor.asin.age.condition.abundance.txt')
CO <- read.delim('CD_Kong/colon/propeller.asin.txt', sep=' ')
TI <- read.delim('CD_Kong/TI/propeller.asin.txt', sep=' ')

cell.prop.lst <- list(pSS=pSS[,c(1,6,8)], UC=UC[,c(1,6,8)], SLE=SLE[,c(1,6,8)], CO=CO[,c(1,6,8)], TI=TI[,c(1,6,8)])

celltypes <- unique(unlist(lapply(cell.prop.lst, function(x) rownames(x))))

# Create matrix
cell.prop.mat <- matrix(0, nrow=length(celltypes), ncol=length(cell.prop.lst))
rownames(cell.prop.mat) <- celltypes
colnames(cell.prop.mat) <- names(cell.prop.lst)

# Fill matrix
for (i in 1:length(cell.prop.lst)) {
    cell.prop.mat[match(cell.prop.lst[[i]][,1], rownames(cell.prop.mat)), i] <- cell.prop.lst[[i]][,2]
}
cell.prop.mat <- cell.prop.mat * -1

# Plot heatmap
pdf('celltype_props_heatmap.pdf')
col=colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
Heatmap(cell.prop.mat, name='T-statistic', col=col, column_title='', row_title='',
column_names_rot=0, row_names_gp=gpar(fontsize=8))
dev.off()