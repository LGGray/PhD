library(edgeR)
library(Seurat)
library(MAST)
library(qvalue)
library(tidyverse)
library(ggplot2)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(reshape2)
library(clusterProfiler)
library(fgsea)

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/Seurat2PB.R')


if(dir.exists('differential.expression/edgeR') != TRUE){dir.create('differential.expression/edgeR')}

# Read in file from command line
pbmc <- readRDS(commandArgs(trailingOnly = TRUE)[1])

# pbmc$age <- as.numeric(gsub('-year-old human stage', '', pbmc$development_stage))

for (cell in levels(pbmc)){
  # Select cell type
  print(cell)
  # subset object by cell type
  pbmc.cell <- subset(pbmc, cellTypist == cell)

  if (length(unique(pbmc.cell$condition)) < 2){
    next
  }

  # Perform Pseudobulking as per edgeR
  y <- Seurat2PB(pbmc.cell, sample='individual', cluster='condition', assay='RNA')

  # Keep samples with library size greater than 1st quartile
  # keep.samples <- y$samples$lib.size > summary(y$samples$lib.size)[2]
  # print('Number of samples dropped and retained')
  # print(table(keep.samples))
  # y <- y[, keep.samples]

  # Remove genes with low expression
  # keep.genes <- filterByExpr(y, group=y$samples$cluster)
  # print('Number of genes dropped and retained')
  # print(table(keep.genes))
  # y <- y[keep.genes, , keep=FALSE]

  # TMM normalisation
  y <- calcNormFactors(y)

  # Reorder cluster i.e condition so disease is the reference
  y$samples$cluster <- factor(y$samples$cluster, levels = c("disease", "control"))
  design <- model.matrix(~ 0 + cluster, data = y$samples)

  # Estimate gene wise dispersion
  y <- estimateDisp(y, design, robust=TRUE)

  # Fit to distribution
  fit <- glmQLFit(y, design, robust=TRUE)

  # Explicitely set up contrasts with disease as reference
  contrastMatrix <- makeContrasts(DiseaseVsControl = clusterdisease - clustercontrol, levels = design)
  
  # Perform QLFTest
  qlf <- glmQLFTest(fit, contrast = contrastMatrix)
  print(summary(decideTests(qlf)))

  # Extract results
  res = topTags(qlf, n = Inf)[[1]]

  # Save to file
  cell = gsub("/|-| ", "_", cell)
  write.table(res, paste0("differential.expression/edgeR/", cell, ".txt"),
              row.names=F, sep="\t", quote = F)
}

print("Done with edgeR-QLF")

# for (cell in levels(pbmc)){
#   # Select cell type
#   print(cell)
#   # subset object by cell type
#   pbmc.cell <- subset(pbmc, cellTypist == cell)

#   # Check if both conditions are present
#   if (length(unique(pbmc.cell$condition)) < 2){
#     next
#   }

#   # Remove genes with pseudobulked expression in less than 5% of individuals
#   expr <- AverageExpression(pbmc.cell, group.by='individual', slot='data')$RNA
#   keep <- apply(expr, 1, function(x) sum(x > 0) > ncol(expr) * 0.05)
#   expr <- expr[keep,]

#   targets = unique(data.frame(group = pbmc.cell$condition,
#                       individual = pbmc.cell$individual))

#   expr.melt <- melt(expr)
#   expr.melt$condition <- targets$group[match(expr.melt$Var2, targets$individual)]
#   expr.melt$condition <- factor(expr.melt$condition, levels = c('disease', 'control'))

#   colnames(expr.melt) <- c('gene', 'individual', 'expression', 'condition')
#   result.list <- lapply(split(expr.melt, expr.melt$gene), function(x){
#     result <- wilcox.test(expression ~ condition, data=x)
#     log2FC <- log2(mean(x[x$condition == 'disease', 'expression']) / mean(x[x$condition == 'control', 'expression']))
#     data.frame(gene=x$gene[[1]], p.value=result$p.value, W=result$statistic, log2FC=log2FC, row.names=NULL)
#   })
#   result.df <- do.call(rbind, result.list)
#   result.df$FDR <- qvalue(p = result.df$p.value)$qvalues

#   result.df$FDR <- p.adjust(result.df$p.value, method='BH')

#   cell = gsub("/|-| ", "_", cell)
#   write.table(result.df, paste0("differential.expression/wilcox/", cell, ".txt"),
#               row.names=F, sep="\t", quote = F)
# }

# ### Analysis of results ###

# source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
# source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')
# load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')

# deg <- deg.list('differential.expression/edgeR', logfc=0.1)

# deg <- deg.list('differential.expression/edgeR', logfc=0.1)
# deg.cellCount <- deg.list('differential.expression/edgeR_cellCount', logfc=0.1)
# deg.age <- deg.list('differential.expression/edgeR_age', logfc=0.1)

# df <- data.frame(celltype=names(deg), deg=unlist(lapply(deg, function(x) nrow(x))),
# deg.age=unlist(lapply(deg.age, function(x) nrow(x))), 
# deg.cellCount=unlist(lapply(deg.cellCount, function(x) nrow(x))))

# # calculate concordance between methods
# df$overlap.age <- unlist(lapply(names(deg), function(x) length(intersect(deg[[x]]$gene, deg.age[[x]]$gene))))
# df$perc.overlap.age <- df$overlap.age / df$deg

# df$overlap.cellCount <- unlist(lapply(names(deg), function(x) length(intersect(deg[[x]]$gene, deg.cellCount[[x]]$gene))))
# df$perc.overlap.cellCount <- df$overlap.cellCount / df$deg

# # round the percentages to 2 decimal places
# df$perc.overlap.age <- round(df$perc.overlap.age * 100, 2)
# df$perc.overlap.cellCount <- round(df$perc.overlap.cellCount * 100, 2)

# df <- df[,c(1,2,3,5,6,4,7,8)]

# df$celltype <- replace.names(gsub('_', '.', df$celltype))

# # Read in propeller results
# propeller <- read.delim('propeller.asin.condition.abundance.txt')
# propeller <- propeller[match(df$celltype, propeller$BaselineProp.clusters),]
# df$Tstatistic <- propeller$Tstatistic * -1
# df$FDR <- propeller$FDR
# write.table(df, 'differential.expression/comparing.nDEG.csv', sep=',', quote=F, row.names=F)

# deg.chrX <- lapply(deg, function(x) subset(x, gene %in% rownames(chrX)))

# ### Heatmap of DEG across celltypes ###
# genes <- unique(unlist(lapply(deg, function(x) x$gene)))
# plot.matrix <- matrix(0, nrow=length(genes), ncol=length(deg))
# rownames(plot.matrix) <- genes
# colnames(plot.matrix) <- replace.names(gsub('_', '.', names(deg)))
# # Match genes to rownames
# for (i in 1:length(deg)){
#   plot.matrix[match(deg[[i]]$gene, genes),i] <- deg[[i]]$logFC.disease_vs_control
# }
# pdf('APR/DEG.heatmap.pdf')
# Heatmap(plot.matrix, clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
# clustering_method_rows = "complete", clustering_method_columns = "complete",
# col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), name='logFC',
# column_title = "Differentially expressed genes", column_title_side = "bottom",
# column_names_rot = 45, column_names_side = "top", column_dend_side = "bottom", show_row_names = FALSE,
# column_names_gp = gpar(fontsize = 9))
# dev.off()

# ### UpSet plot of DEG ###
# deg.lst <- lapply(deg, function(x) x$gene)
# deg.mtx <- fromList(deg.lst)
# rownames(deg.mtx) <- unique(unlist(deg.lst))
# colnames(deg.mtx) <- replace.names(gsub('_', '.', names(deg)))
# pdf('APR/DEG.upset.pdf', width=12, height=10)
# upset(deg.mtx, order.by = "freq", nsets = length(deg.lst), nintersects=NA, 
#       main.bar.color = "black", sets.bar.color = "black", 
#       matrix.color = "black", shade.color = "black")
# dev.off()

# merged <- merge(deg$HSC_MPP, deg$Plasma_cells, by='gene')
# cor.test(merged$logFC.disease_vs_control.x, merged$logFC.disease_vs_control.y, method='spearman')

# ### Heatmap of chrX genes across celltypes ###
# genes.chrX <- unique(unlist(lapply(deg.chrX, function(x) x$gene)))
# plot.matrix.chrX <- matrix(0, nrow=length(genes.chrX), ncol=length(deg.chrX))
# rownames(plot.matrix.chrX ) <- genes.chrX
# colnames(plot.matrix.chrX ) <- replace.names(gsub('_', '.', names(deg.chrX)))
# # Match genes to rownames
# for (i in 1:length(deg.chrX)){
#   if (nrow(deg.chrX[[i]]) > 0) {
#     plot.matrix.chrX[match(deg.chrX[[i]]$gene, genes.chrX), i] <- deg.chrX[[i]]$logFC.disease_vs_control
#   }
# }
# pdf('APR/DEG.chrX.heatmap.pdf')
# Heatmap(plot.matrix.chrX, clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
# clustering_method_rows = "complete", clustering_method_columns = "complete",
# col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), name='logFC', 
# column_title = "Differentially expressed X chromosome genes", column_title_side = "bottom",
# column_names_rot = 45, column_names_side = "top", column_dend_side = "bottom", show_row_names = FALSE,
# column_names_gp = gpar(fontsize = 9))
# dev.off()

# ### UpSet plot of chrX DEG ###
# deg.chrX.lst <- lapply(deg.chrX, function(x) x$gene)
# deg.chrX.mtx <- fromList(deg.chrX.lst)
# rownames(deg.chrX.mtx) <- unique(unlist(deg.chrX.lst))
# colnames(deg.chrX.mtx) <- replace.names(gsub('_', '.', names(deg.chrX)))
# pdf('APR/DEG.chrX.upset.pdf', width=10, height=10)
# upset(deg.chrX.mtx, order.by = "freq", nsets = length(deg.chrX.lst), nintersects=NA, point.size = 2, line.size = 1.5, 
#       main.bar.color = "black", sets.bar.color = "black", text.scale = 1.5, 
#       matrix.color = "black", shade.color = "black")
# dev.off()

# rownames(deg.chrX.mtx)[order(rowSums(deg.chrX.mtx), decreasing = TRUE)[1:10]]

# ### Calculating enrichment ###
# source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/fishers.test.degs.R')
# edgeR <- deg.list('differential.expression/edgeR_cellCount', filter=F)
# names(edgeR) <- replace.names(gsub('_', '.', names(edgeR)))

# ### Fishers test chrX ###
# # all
# fishers.all.chrX <- lapply(edgeR, function(x) fisher.test.edgeR(x, rownames(chrX), 0.1, direction='none'))
# fishers.all.chrX[sapply(fishers.all.chrX, function(x) x$p.value < 0.05)]

# # up
# fishers.up.chrX <- lapply(edgeR, function(x) fisher.test.edgeR(x, rownames(chrX), 0.1, direction='up'))
# fishers.up.chrX[sapply(fishers.up.chrX, function(x) x$p.value < 0.05)]

# # down
# fishers.down.chrX <- lapply(edgeR, function(x) fisher.test.edgeR(x, rownames(chrX), 0.1, direction='down'))
# fishers.down.chrX[sapply(fishers.down.chrX, function(x) x$p.value < 0.05)]

# ### Fishers test DMG ###
# # all
# load('/directflow/SCCGGroupShare/projects/lacgra/datasets/DMG.Rdata')
# fishers.all.DMG <- lapply(edgeR, function(x) fisher.test.edgeR(x, DMG, 0.1, direction='none'))
# fishers.all.DMG[sapply(fishers.all.DMG, function(x) x$p.value < 0.05)]

# # up
# fishers.up.DMG <- lapply(edgeR, function(x) fisher.test.edgeR(x, DMG, 0.1, direction='up'))
# fishers.up.DMG[sapply(fishers.up.DMG, function(x) x$p.value < 0.05)]

# # down
# fishers.down.DMG <- lapply(edgeR, function(x) fisher.test.edgeR(x, DMG, 0.1, direction='down'))
# fishers.down.DMG[sapply(fishers.down.DMG, function(x) x$p.value < 0.05)]

# # Plot heatmap of DMG across celltypes
# deg.DMG <- lapply(deg, function(x) subset(x, gene %in% DMG))
# ### Heatmap of chrX genes across celltypes ###
# genes.DMG <- unique(unlist(lapply(deg.DMG, function(x) x$gene)))
# plot.matrix.DMG<- matrix(0, nrow=length(genes.DMG), ncol=length(deg.DMG))
# rownames(plot.matrix.DMG) <- genes.DMG
# colnames(plot.matrix.DMG) <- replace.names(gsub('_', '.', names(deg.DMG)))
# # Match genes to rownames
# for (i in 1:length(deg.DMG)){
#   if (nrow(deg.DMG[[i]]) > 0) {
#     plot.matrix.DMG[match(deg.DMG[[i]]$gene, genes.DMG), i] <- deg.DMG[[i]]$logFC.disease_vs_control
#   }
# }
# pdf('APR/DEG.DMG.heatmap.pdf')
# Heatmap(plot.matrix.DMG, clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
# clustering_method_rows = "complete", clustering_method_columns = "complete",
# col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), name='logFC', 
# column_title = "Differentially expressed differentially methylated genes", column_title_side = "bottom",
# column_names_rot = 45, column_names_side = "top", column_dend_side = "bottom", show_row_names = TRUE,
# column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9))
# dev.off()

# ### Fishers test SLE DisGeNet ###
# SLE <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/SLE.tsv')
# # all
# fishers.all.SLE <- lapply(edgeR, function(x) fisher.test.edgeR(x, SLE$Gene, 0.1, direction='none'))
# names(fishers.all.SLE[sapply(fishers.all.SLE, function(x) x$p.value < 0.05)])

# # up
# fishers.up.SLE <- lapply(edgeR, function(x) fisher.test.edgeR(x, SLE$Gene, 0.1, direction='up'))
# fishers.up.SLE[sapply(fishers.up.SLE, function(x) x$p.value < 0.05)]

# # down
# fishers.down.SLE <- lapply(edgeR, function(x) fisher.test.edgeR(x, SLE$Gene, 0.1, direction='down'))
# fishers.down.SLE[sapply(fishers.down.SLE, function(x) x$p.value < 0.05)]

# # Plot heatmap of SLE across celltypes
# deg.SLE <- lapply(deg, function(x) subset(x, gene %in% SLE$Gene))
# ### Heatmap of chrX genes across celltypes ###
# genes.SLE <- unique(unlist(lapply(deg.SLE, function(x) x$gene)))
# plot.matrix.SLE<- matrix(0, nrow=length(genes.SLE), ncol=length(deg.SLE))
# rownames(plot.matrix.SLE) <- genes.SLE
# colnames(plot.matrix.SLE) <- replace.names(gsub('_', '.', names(deg.SLE)))
# # Match genes to rownames
# for (i in 1:length(deg.SLE)){
#   if (nrow(deg.SLE[[i]]) > 0) {
#     plot.matrix.SLE[match(deg.SLE[[i]]$gene, genes.SLE), i] <- deg.SLE[[i]]$logFC.disease_vs_control
#   }
# }
# pdf('APR/DEG.SLE.heatmap.pdf')
# Heatmap(plot.matrix.SLE, clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
# clustering_method_rows = "complete", clustering_method_columns = "complete",
# col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), name='logFC',
# column_title = "Differentially expressed SLE genes", column_title_side = "bottom",
# column_names_rot = 45, column_names_side = "top", column_dend_side = "bottom", show_row_names = FALSE)
# dev.off()

# ### UpSet plot of SLE DEG ###
# deg.SLE.lst <- lapply(deg.SLE, function(x) x$gene)
# deg.SLE.mtx <- fromList(deg.SLE.lst)
# rownames(deg.SLE.mtx) <- unique(unlist(deg.SLE.lst))
# colnames(deg.SLE.mtx) <- replace.names(gsub('_', '.', names(deg.SLE)))

# pdf('APR/DEG.SLE.upset.pdf', width=10, height=10)
# upset(deg.SLE.mtx, order.by = "freq", nsets = length(deg.SLE.lst), nintersects=NA, point.size = 2, line.size = 1.5, 
#       main.bar.color = "black", sets.bar.color = "black", text.scale = 1.5, 
#       matrix.color = "black", shade.color = "black")
# dev.off()


# ### fgsea of hallmark, GO, KEGG and Reactome ###
# library(msigdbr)
# hallmark = msigdbr(species = "human", category = "H")
# hallmark_list = split(x = hallmark$gene_symbol, f = hallmark$gs_name)
# GOBP = msigdbr(species = "human", category = "C5", subcategory = "BP")
# GOBP_list = split(x = GOBP$gene_symbol, f = GOBP$gs_name)
# KEGG = msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
# KEGG_list = split(x = KEGG$gene_symbol, f = KEGG$gs_name)
# Reactome = msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME")
# Reactome_list = split(x = Reactome$gene_symbol, f = Reactome$gs_name)
# pathways <- list(hallmark_list, GOBP_list, KEGG_list, Reactome_list)

# for(cell in names(edgeR)){
#   print(cell)
#   ranked_gene_list <- edgeR[[cell]]$logFC.disease_vs_control
#   names(ranked_gene_list) <- edgeR[[cell]]$gene
#   # Rank genes based on the absolute value of logFC
#   ranked_gene_list <- ranked_gene_list[order(abs(ranked_gene_list), decreasing = TRUE)]
  
#   fgsea_res <- lapply(pathways, function(x){
#     res <- fgsea(pathways = x, stats = ranked_gene_list, minSize  = 15, maxSize  = 500)
#     collapsedPathways <- collapsePathways(res[order(pval)][padj < 0.01], x, ranked_gene_list)
#     res[pathway %in% collapsedPathways$mainPathways][order(-NES),]
#   })
#   names(fgsea_res) <- c('hallmark', 'GOBP', 'KEGG', 'Reactome')
  
#   fgsea_final <- bind_rows(fgsea_res, .id='geneSet')
#   data.table::fwrite(fgsea_final, file=paste0('fgsea/', gsub("/|-| ", "_", cell), '.txt'), sep="\t", sep2=c("", " ", ""))

#   # Order the pathways by padj
#   fgsea_final <- fgsea_final[order(fgsea_final$padj), ]
#   fgsea_final$pathway <- factor(fgsea_final$pathway, levels = unique(fgsea_final$pathway))
#   fgsea_final$geneSet <- factor(fgsea_final$geneSet, levels = c('hallmark', 'GOBP', 'KEGG', 'Reactome'))
#   fgsea_final$chrX <- unlist(lapply(fgsea_final$leadingEdge, function(x) length(intersect(x, rownames(chrX)))))

#   # Create the plot
#   pdf(paste0('fgsea/', gsub("/|-| ", "_", cell), '.pdf'), width=20, height=15)
#   p <- ggplot(fgsea_final, aes(x=pathway, y=-log10(padj), color=geneSet, size=size, alpha=chrX)) +
#   geom_point() +
#   coord_flip() +
#   ggtitle(cell) +
#   theme(axis.text.x = element_text(size=5)) +
#   scale_size_continuous(name = "# genes in pathway") +
#   scale_alpha_continuous(name = "# chrX genes in pathway", range = c(0.5, 1)) +
#   scale_color_discrete(name = "geneSet")
#   print(p)
#   dev.off()