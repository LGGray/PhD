### AUCell of XCI escape genes ###
library(AUCell)
library(Seurat)
library(GSEABase)

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')

pbmc <- readRDS('pbmc.female.control-managed.RDS')

# save all degs to file for input into interferomedb
degs <- deg.list('differential.expression/edgeR', logfc=0.1)
gene.list <- unique(unlist(lapply(degs, function(x){x$gene})))
write.csv(gene.list, file='all.degs.csv', row.names=F)

# Read in interferome results
interferome <- unique(read.delim('interferomedb_SLE.txt', skip=19)$Gene.Name)

# Run AUCell
if(dir.exists('AUCell')==FALSE){
    dir.create('AUCell')
}
geneSets <- GeneSet(interferome, setName="ISG")

exprMatrix <- GetAssayData(pbmc, assay="RNA", slot="data")

cells_AUC <- AUCell_run(exprMatrix, geneSets, BPPARAM=BiocParallel::MulticoreParam(4))

pdf('AUCell/AUCell_thresholds.pdf')
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
dev.off()

pbmc@meta.data[rownames(pbmc@meta.data) %in% cells_assignment$ISG$assignment,c('condition','cellTypist')]

auc_matrix <- getAUC(cells_AUC)[,cells_assignment$ISG$assignment]

pbmc <- AddMetaData(pbmc, auc_matrix, col.name='AUCell', assay='RNA')


meta <- pbmc@meta.data

wilcox.results <- list()
plot.data.list <- list()
for(cell in levels(pbmc)){
    control_list <- lapply(unique(subset(meta, condition=='control')$individual), function(x) {
    cells <- rownames(subset(meta, individual==x & cellTypist==cell & condition=='control'))
    auc <- mean(auc_matrix[1,cells])
    })
    control <- unlist(control_list)
    disease_list <- lapply(unique(subset(meta, condition=='disease')$individual), function(x) {
        cells <- rownames(subset(meta, individual==x & cellTypist==cell & condition=='disease'))
        auc <- mean(auc_matrix[1,cells])
    })
    disease <- unlist(disease_list)

    plot.data <- rbind(data.frame(condition='control', AUC=control), data.frame(condition='disease', AUC=disease))
    plot.data$condition <- factor(plot.data$condition, levels=c('disease', 'control'))
    plot.data.list[[cell]] <- plot.data

    tmp <- wilcox.test(AUC ~ condition, data=plot.data, alternative='greater')
    wilcox.results[[cell]] <- tmp$p.value
}

wilcox.result.combined <- data.frame(celltype=names(wilcox.results), pvalue=unlist(wilcox.results))
wilcox.result.combined$FDR <- p.adjust(wilcox.result.combined$pvalue, method='fdr')
wilcox.result.combined <- wilcox.result.combined[order(wilcox.result.combined$FDR, decreasing=TRUE),]
wilcox.result.combined$celltype <- factor(wilcox.result.combined$celltype, levels=wilcox.result.combined$celltype)
write.table(wilcox.result.combined, file='AUCell/wilcox.result.combined.txt', sep='\t', row.names=F)

pdf('AUCell/AUCell_escape.pdf')
ggplot(wilcox.result.combined, aes(x=celltype, y=-log10(FDR), fill=-log10(FDR))) +
geom_bar(stat='identity') +
scale_fill_gradient(low='blue', high='red') +
geom_hline(yintercept=-log10(0.05), linetype='dashed') +
theme(legend.position='none') +
xlab('') +
coord_flip()
dev.off()


pdf('AUCell/AUCell_escape.violin.pdf')
ggplot(plot.data.list[['Tem/Temra cytotoxic T cells']], aes(x=condition, y=AUC, fill=condition)) + 
geom_violin() +
scale_fill_manual(values=c('disease'='red', 'control'='blue')) +
geom_boxplot(width=0.1)
dev.off()

save(auc_matrix, file='AUCell_katsir.RData')