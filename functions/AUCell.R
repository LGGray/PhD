### AUCell of XCI escape genes ###
library(AUCell)
library(Seurat)

pbmc <- readRDS('../pbmc.female.control-managed.RDS')
disease <- 'SLE'

# load in gene sets
escape <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/Katsir.escape.txt')
escape <- escape$Gene.Symbol
chrX <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX_biomaRt.txt')
chrX <- subset(chrX, Gene.name != '')
chrX <- unique(chrX$Gene.name)
hallmark <- clusterProfiler::read.gmt('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/h.all.v7.5.1.symbols.gmt')
estrogen <- unique(subset(hallmark, term %in% c('HALLMARK_ESTROGEN_RESPONSE_EARLY','HALLMARK_ESTROGEN_RESPONSE_LATE'))$gene)
androgen <- unique(subset(hallmark, term %in% c('HALLMARK_ANDROGEN_RESPONSE'))$gene)
disgene <- read.delim(paste0('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/', disease, '.tsv'))$Gene
GWAS <- read.delim(paste0('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/', disease, '_GWAS.tsv'))
GWAS <- unique(unlist(lapply(GWAS$Gene, function(x) unlist(strsplit(x, ';')))))

geneSets <- list(escape=escape, chrX=chrX, estrogen=estrogen, androgen=androgen, disgene=disgene, GWAS=GWAS)

cell <- 'CD16+ NK cells'
for(cell in levels(pbmc)){
    pbmc.subset <- subset(pbmc, cellTypist==cell)

    exprMatrix <- GetAssayData(pbmc.subset, assay="RNA", slot="data")

    # Rank genes in each cell
    cells_rankings <- AUCell_buildRankings(exprMatrix)
    saveRDS(cells_rankings, paste0(gsub(' ', '_', cell),'_rankings.RDS'))

    cells_AUC <- AUCell_run(exprMatrix, geneSets, BPPARAM=BiocParallel::MulticoreParam(4))
    save(cells_AUC, file='cells_AUC.RData')

    pdf('AUCell_thresholds.pdf')
    par(mfrow=c(3,2))
    cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
    dev.off()

    auc_matrix <- getAUC(cells_AUC)
    save(auc_matrix, file='auc_matrix.RData')
}

cells_assignment$escape$aucThr$thresholds


# meta <- pbmc@meta.data

# wilcox.results <- list()
# plot.data.list <- list()
# for(cell in levels(pbmc)){
#     control_list <- lapply(unique(subset(meta, condition=='control')$individual), function(x) {
#     cells <- rownames(subset(meta, individual==x & cellTypist==cell & condition=='control'))
#     auc <- mean(auc_matrix[1,cells])
#     })
#     control <- unlist(control_list)
#     disease_list <- lapply(unique(subset(meta, condition=='disease')$individual), function(x) {
#         cells <- rownames(subset(meta, individual==x & cellTypist==cell & condition=='disease'))
#         auc <- mean(auc_matrix[1,cells])
#     })
#     disease <- unlist(disease_list)

#     plot.data <- rbind(data.frame(condition='control', AUC=control), data.frame(condition='disease', AUC=disease))
#     plot.data$condition <- factor(plot.data$condition, levels=c('disease', 'control'))
#     plot.data.list[[cell]] <- plot.data

#     tmp <- wilcox.test(AUC ~ condition, data=plot.data, alternative='greater')
#     wilcox.results[[cell]] <- tmp$p.value
# }

# wilcox.result.combined <- data.frame(celltype=names(wilcox.results), pvalue=unlist(wilcox.results))
# wilcox.result.combined$FDR <- p.adjust(wilcox.result.combined$pvalue, method='fdr')
# wilcox.result.combined <- wilcox.result.combined[order(wilcox.result.combined$FDR, decreasing=TRUE),]
# wilcox.result.combined$celltype <- factor(wilcox.result.combined$celltype, levels=wilcox.result.combined$celltype)
# write.table(wilcox.result.combined, file='AUCell/wilcox.result.combined.txt', sep='\t', row.names=F)

# pdf('AUCell/AUCell_escape.pdf')
# ggplot(wilcox.result.combined, aes(x=celltype, y=-log10(FDR), fill=-log10(FDR))) +
# geom_bar(stat='identity') +
# scale_fill_gradient(low='blue', high='red') +
# geom_hline(yintercept=-log10(0.05), linetype='dashed') +
# theme(legend.position='none') +
# xlab('') +
# coord_flip()
# dev.off()


# pdf('AUCell/AUCell_escape.violin.pdf')
# ggplot(plot.data.list[['Tem/Temra cytotoxic T cells']], aes(x=condition, y=AUC, fill=condition)) + 
# geom_violin() +
# scale_fill_manual(values=c('disease'='red', 'control'='blue')) +
# geom_boxplot(width=0.1)
# dev.off()

# save(auc_matrix, file='AUCell_katsir.RData')