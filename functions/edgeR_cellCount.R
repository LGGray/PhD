library(edgeR)
library(Seurat)
library(qvalue)
library(tidyverse)
library(ggplot2)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(reshape2)
library(clusterProfiler)

if(dir.exists('differential.expression') != TRUE){dir.create('differential.expression')}
if(dir.exists('differential.expression/edgeR') != TRUE){dir.create('differential.expression/edgeR')}
if(dir.exists('differential.expression/MAST') != TRUE){dir.create('differential.expression/MAST')}

# Read in file from command line
pbmc <- readRDS(commandArgs(trailingOnly = TRUE)[1])
assay <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

# unique(paste(pbmc$condition, pbmc$individual, pbmc$Type))
# pbmc <- subset(pbmc, Type %in% c('Heal', 'NonI'))

pbmc$age <- as.numeric(gsub('-year-old human stage', '', pbmc$development_stage))

for (cell in levels(pbmc)){
  # Select cell type
  print(cell)
  # subset object by cell type
  pbmc.cell <- subset(pbmc, cellTypist == cell)
  # check if there are enough cells and skip if not
  if(nrow(pbmc.cell) < 10){
    print("Not enough cells")
    next
  }
  # check if there are enough cell in both conditions and skip if not
  if(length(unique(pbmc.cell$condition)) != 2){
    print("Not enough conditions")
    next
  } else {
    table(pbmc.cell$condition, pbmc.cell$cellTypist)
  }
  # Keep genes with expression in 5% of cells
  keep <- rowSums(pbmc.cell@assays$RNA@counts > 0) > ncol(pbmc.cell) * 0.05
  features <- names(keep[keep == T])
  pbmc.cell <- subset(pbmc.cell, features=features)

  # Calculate cellcount
  cellCount <- as.data.frame.matrix(table(pbmc.cell$individual, pbmc.cell$cellTypist))

  # Psudobulking by summing counts
  expr <- AggregateExpression(pbmc.cell, group.by='individual', slot='counts')[[assay]]
  expr <- expr[,(colSums(expr) > 0)]

  # edgeR-QLFTest
  targets = unique(data.frame(condition = pbmc.cell$condition,
                      individual = pbmc.cell$individual,
                      age = pbmc.cell$age))
  targets$cellCount <- cellCount[,cell]
  targets <- targets[match(colnames(expr), targets$individual),]

  # targets$SV1 <- tapply(pbmc.cell$SV1, pbmc.cell$individual, sum)
  # targets$SV2 <- tapply(pbmc.cell$SV2, pbmc.cell$individual, sum)
  rownames(targets) <- targets$individual
  design <- model.matrix(~0 + cellCount + age + condition, data=targets)
  y = DGEList(counts = expr, group = targets$condition)
  contrasts <- makeContrasts(disease_vs_control = conditiondisease - conditioncontrol,
                            cell_type_effect = cellCount,
                            levels = design) 
  # Disease group as reference
  y = DGEList(counts = expr, group = targets$condition)
  y = estimateGLMRobustDisp(y, design, trend.method = 'auto')
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, contrast=contrasts)
  print(summary(decideTests(qlf)))
    res = topTags(qlf, n = Inf) %>%
      as.data.frame() %>%
      rownames_to_column('gene')
    tryCatch({
      res$FDR <- qvalue(p = res$PValue)$qvalues
    }, error = function(e) {
      print(paste0("Error in qvalue(): ", e$message))
    })
  cell = gsub("/|-| ", "_", cell)
  write.table(res, paste0("differential.expression/edgeR/", cell, ".txt"),
              row.names=F, sep="\t", quote = F)
}

print("Done with edgeR-QLF")

### Analysis of results ###

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/replace.names.R')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')

edgeR <- deg.list('differential.expression/edgeR', filter=F)
deg <- deg.list('differential.expression/edgeR', logfc=0.5)
# names(deg) <- c(
#   "CD16+ NK cells", "Classical monocytes", "DC1", "DC2", "MAIT cells", "Mast cells", "Memory B cells", 
#   "Naive B cells", "NK cells", "Non-classical monocytes", "pDC", "Plasma cells", "Plasmablasts", "Regulatory T cells", 
#   "Tcm/Naive cytotoxic T cells", "Tcm/Naive helper T cells", "Tem/Effector helper T cells", "Tem/Temra cytotoxic T cells", 
#   "Tem/Trm cytotoxic T cells"
# )

# names(deg) <- c(
#   "CD16+ NK cells", "Classical monocytes", "CRTAM+ gamma delta T cells", "DC2", "gamma delta T cells", "ILC3", "MAIT cells", "Memory B cells", 
#   "Naive B cells", "Non-classical monocytes", "Plasma cells", "Regulatory T cells", "Tcm/Naive cytotoxic T cells", "Tcm/Naive helper T cells", 
#   "Tem/Effector helper T cells", "Tem/Temra_cytotoxic T cells", "Tem/Trm cytotoxic T cells", "Trm/cytotoxic T cells", "Type 1 helper T cells"
# )

cell_types <- c("CD16- NK_cells", "CD16+ NK_cells", "Classical monocytes", "Cycling T cells", "DC1", "DC2", "Erythrophagocytic macrophages", "Follicular_helper_T_cells",
        "gamma delta T cells", "Germinal center B cells", "ILC", "Intestinal macrophages", "Macrophages", "MAIT cells", "Mast cells", "Memory B cells", 
        "Migratory DCs", "Myelocytes", "Naive B cells", "Plasma cells", "Proliferative germinal center B cells","Regulatory T cells", "Tcm/Naive helper T cells", 
        "Tem/Effector helper T cells", "Tem/Trm_cytotoxic T cells", "Trm/cytotoxic T cells", "Type/17 helper T cells")
names(deg) <- cell_types

deg.chrX <- lapply(deg, function(x) subset(x, gene %in% rownames(chrX)))

# Heatmap of all genes across celltypes
genes <- unique(unlist(lapply(edgeR, function(x) x$gene)))
plot.matrix <- matrix(0, nrow=length(genes), ncol=length(edgeR))
rownames(plot.matrix) <- genes
colnames(plot.matrix) <- replace.names(gsub('_', '.', names(edgeR)))
# Match genes to rownames
for (i in 1:length(edgeR)){
  plot.matrix[match(edgeR[[i]]$gene, genes),i] <- edgeR[[i]]$logFC.disease_vs_control
}
pdf('APR/all.genes.heatmap.pdf')
Heatmap(plot.matrix, clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
clustering_method_rows = "complete", clustering_method_columns = "complete",
col=colorRamp2(c(-5, 0, 5), c("blue", "white", "red")), name='logFC',
column_title = "All genes", column_title_side = "bottom",
column_names_rot = 45, column_names_side = "top", column_dend_side = "bottom", show_row_names = FALSE)
dev.off()


# Heatmap of DEG across celltypes
genes <- unique(unlist(lapply(deg, function(x) x$gene)))
plot.matrix <- matrix(0, nrow=length(genes), ncol=length(deg))
rownames(plot.matrix) <- genes
colnames(plot.matrix) <- replace.names(gsub('_', '.', names(deg)))
# Match genes to rownames
for (i in 1:length(deg)){
  plot.matrix[match(deg[[i]]$gene, genes),i] <- deg[[i]]$logFC.disease_vs_control
}
pdf('APR/DEG.heatmap.pdf')
Heatmap(plot.matrix, clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
clustering_method_rows = "complete", clustering_method_columns = "complete",
col=colorRamp2(c(-5, 0, 5), c("blue", "white", "red")), name='logFC',
column_title = "Differentially expressed genes", column_title_side = "bottom",
column_names_rot = 45, column_names_side = "top", column_dend_side = "bottom", show_row_names = FALSE,
column_names_gp = gpar(fontsize = 9))
dev.off()

# Correlation of DEG
plot.matrix.cor <- cor(plot.matrix, method='spearman')
pdf('APR/DEG.heatmap.cor.pdf')
Heatmap(plot.matrix.cor, clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
clustering_method_rows = "complete", clustering_method_columns = "complete",
col=colorRamp2(c(0, 0.5, 1), c("blue","white","red")), name='Rho',
column_title = "Correlation of DEG between Cell Types",  column_title_side = "bottom",
column_names_rot = 45, column_names_side = "top", column_dend_side = "bottom", show_row_names = FALSE,
column_names_gp = gpar(fontsize = 9))
dev.off()

# Heatmap of chrX genes across celltypes
genes.chrX <- unique(unlist(lapply(deg.chrX, function(x) x$gene)))
plot.matrix.chrX <- matrix(0, nrow=length(genes.chrX), ncol=length(deg.chrX))
rownames(plot.matrix.chrX ) <- genes.chrX
colnames(plot.matrix.chrX ) <- replace.names(gsub('_', '.', names(deg.chrX)))
# Match genes to rownames
for (i in 1:length(deg.chrX)){
  if (nrow(deg.chrX[[i]]) > 0) {
    plot.matrix.chrX[match(deg.chrX[[i]]$gene, genes.chrX), i] <- deg.chrX[[i]]$logFC.disease_vs_control
  }
}

pdf('APR/DEG.chrX.heatmap.pdf')
Heatmap(plot.matrix.chrX, clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
clustering_method_rows = "complete", clustering_method_columns = "complete",
col=colorRamp2(c(-5, 0, 5), c("blue", "white", "red")), name='logFC', 
column_title = "Differentially expressed X chromosome genes", column_title_side = "bottom",
column_names_rot = 45, column_names_side = "top", column_dend_side = "bottom", show_row_names = FALSE,
column_names_gp = gpar(fontsize = 9))
dev.off()

# Correlation of chrX DEG
# Remove cells with variance == 0
plot.matrix.chrX <- plot.matrix.chrX[,apply(plot.matrix.chrX, 2, var) > 0]
plot.matrix.chrX.cor <- cor(plot.matrix.chrX, method='spearman')
pdf('APR/DEG.chrX.heatmap.cor.pdf')
Heatmap(plot.matrix.chrX.cor, clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
clustering_method_rows = "complete", clustering_method_columns = "complete",
col=colorRamp2(c(0, 0.5, 1), c("blue","white","red")), name='Rho',
column_title = "Correlation of DEG chrX",  column_title_side = "bottom",
column_names_rot = 45, column_names_side = "top", column_dend_side = "bottom", show_row_names = FALSE,
column_names_gp = gpar(fontsize = 9))
dev.off()

### UpSet plot of DEG ###
deg.lst <- lapply(deg, function(x) x$gene)
deg.mtx <- fromList(deg.lst)
colnames(deg.mtx) <- replace.names(gsub('_', '.', names(deg)))
pdf('APR/DEG.upset.pdf')
upset(deg.mtx, order.by = "freq", nsets = length(deg.lst), nintersects=NA, point.size = 2, line.size = 1.5, 
      main.bar.color = "black", sets.bar.color = "black", text.scale = 1.5, 
      matrix.color = "black", shade.color = "black")
dev.off()

### UpSet plot of chrX DEG ###
deg.chrX.lst <- lapply(deg.chrX, function(x) x$gene)
deg.chrX.mtx <- fromList(deg.chrX.lst)
colnames(deg.chrX.mtx) <- replace.names(gsub('_', '.', names(deg.chrX)))
pdf('APR/DEG.chrX.upset.pdf', width=10, height=10)
upset(deg.chrX.mtx, order.by = "freq", nsets = length(deg.chrX.lst), nintersects=NA, point.size = 2, line.size = 1.5, 
      main.bar.color = "black", sets.bar.color = "black", text.scale = 1.5, 
      matrix.color = "black", shade.color = "black")
dev.off()




### Calculating enrichment ###

edgeR <- deg.list('differential.expression/edgeR', filter=F)
names(edgeR) <- replace.names(gsub('_', '.', names(edgeR)))

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/chisq.test.degs.R')
chisq.up <- lapply(edgeR, function(x) chisq.test.edgeR(x, rownames(chrX), 0.5, direction='up'))
chisq.up.df <- dplyr::bind_rows(lapply(chisq.up, function(x) data.frame(pvalue=x$p.value, statistic=x$statistic[[1]])), .id='celltype')
chisq.up.df$size <- unlist(lapply(deg, function(x) nrow(subset(x, gene %in% rownames(chrX) & logFC.disease_vs_control > 0.5))))
write.table(chisq.up.df, 'APR/chrX.up.chisq.txt', sep='\t', row.names=F)

chisq.down <- lapply(edgeR, function(x) chisq.test.edgeR(x, rownames(chrX), 0.5, direction='down'))
chisq.down.df <- dplyr::bind_rows(lapply(chisq.down, function(x) data.frame(pvalue=x$p.value, statistic=x$statistic[[1]])), .id='celltype')
chisq.down.df$size <- unlist(lapply(deg, function(x) nrow(subset(x, gene %in% rownames(chrX) & logFC.disease_vs_control < -0.5))))
write.table(chisq.down.df, 'APR/chrX.down.chisq.txt', sep='\t', row.names=F)

# Plot the results
chisq.up.df <- chisq.up.df[chisq.up.df$size > 0,]
pdf('APR/chrX.up.chisq.pdf')
ggplot(chisq.up.df, aes(x=statistic, y=-log10(pvalue))) + 
    geom_point(aes(size=size, fill=log10(pvalue))) +
    geom_text(aes(label=celltype)) +
    ylab("-log10(pvalue)") + xlab("Chi-squared statistic") + 
    ggtitle("Upregulated X chromosome genes") +
    theme_bw() + 
    theme(plot.title = element_text(size = 18, hjust = 0))
dev.off()

# Fishers test - all
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/fishers.test.degs.R')
fishers.all <- lapply(edgeR, function(x) fisher.test.edgeR(x, rownames(chrX), 0.5, direction='none'))

# Fishers test - up
fishers.up <- lapply(edgeR, function(x) fisher.test.edgeR(x, rownames(chrX), 0.5, direction='up'))
names(fishers.up) <- names(edgeR)
fishers.up.df <- dplyr::bind_rows(lapply(fishers.up, function(x) data.frame(pvalue=x$p.value, statistic=x$estimate[[1]])), .id='celltype')
fishers.up.df$size <- unlist(lapply(deg, function(x) nrow(subset(x, gene %in% rownames(chrX) & logFC.disease_vs_control > 0.5))))
fishers.up.df
write.table(fishers.up.df, 'APR/chrX.up.fishers.txt', sep='\t', row.names=F)
# Fishers test - down
fishers.down <- lapply(edgeR, function(x) fisher.test.edgeR(x, rownames(chrX), 0.5, direction='down'))
fishers.down.df <- dplyr::bind_rows(lapply(fishers.down, function(x) data.frame(pvalue=x$p.value, statistic=x$estimate[[1]])), .id='celltype')
fishers.down.df$size <- unlist(lapply(deg, function(x) nrow(subset(x, gene %in% rownames(chrX) & logFC.disease_vs_control < -0.5))))
fishers.down.df


### Enrichment of genes in hallmark pathways ###
# Read in hallmark pathways
hallmark <- read.gmt('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/h.all.v7.5.1.symbols.gmt')

pathway.enrich <- function(results, universe, gene.set){
  a <- length(results[results %in% gene.set])
  b <- length(!(universe[universe %in% gene.set] %in% results))
  c <- length(results[!(results %in% gene.set)])
  d <- length(!(universe[!(universe %in% gene.set)] %in% results))
  fisher.test(matrix(c(a, b, c, d), nrow=2), alternative='greater')$p.value
}

pathway.analysis <- list()
for(cell in names(edgeR)){
  df <- edgeR[[cell]]
  DEG <- subset(df, FDR < 0.05 & logFC.disease_vs_control > 0.5)$gene

  enrich <- enricher(gene = DEG,
    universe = df$gene,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    TERM2GENE = hallmark)
  if(is.null(enrich)){next}
  enrich <- subset(enrich@result, p.adjust < 0.05)
  # Check which pathways include X chromosome genes
  enrich.chrX <- lapply(enrich$geneID, function(x){
    tmp <- unlist(strsplit(x, '/'))
    pvalue <- pathway.enrich(tmp, DEG, rownames(chrX))
    list(genes=tmp[tmp %in% rownames(chrX)], pvalue=pvalue)
  })
  names(enrich.chrX) <- enrich$ID
  if(length(enrich.chrX) == 0){
    enrich.chrX <- list()
    pathway.analysis[[cell]] <- list(all=enrich, chrX=enrich.chrX)
  } else {
  enrich.chrX <- enrich.chrX[sapply(enrich.chrX, function(x) length(x$genes) > 0)]
  pathway.analysis[[cell]] <- list(all=enrich, chrX=enrich.chrX)
  }
}

# Create matrix of results
all.list <- fromList(lapply(pathway.analysis, function(x) x$all$ID))
rownames(all.list) <- unique(unlist(lapply(pathway.analysis, function(x) x$all$ID)))

chrX.list <- fromList(lapply(pathway.analysis, function(x) names(x$chrX)))
unlist(lapply(pathway.analysis, function(x) lapply(x$chrX, function(y) y$pvalue)))
rownames(chrX.list) <- unique(unlist(lapply(pathway.analysis, function(x) names(x$chrX))))

# Heatmap of pathways
pdf('APR/all.up.Hallmark.heatmap.pdf', width=12, height=10)
Heatmap(as.matrix(all.list), name='mat', col=colorRamp2(c(0, 1), c("white", "red")),
row_names_side = "left", cluster_rows = FALSE, show_column_dend = FALSE, column_names_rot = 45, width = unit(12, "cm"),
column_title = "Upregulated pathways")
dev.off()

pdf('APR/chrX.up.Hallmark.heatmap.pdf', width=12, height=10)
Heatmap(as.matrix(chrX.list), name='mat', col=colorRamp2(c(0, 1), c("white", "red")),
row_names_side = "left", cluster_rows = FALSE, show_column_dend = FALSE, column_names_rot = 45, width = unit(10, "cm"),
column_title = "Upregulated chrX pathways")
dev.off()


pdf('APR/all.down.Hallmark.heatmap.pdf', width=12, height=10)
Heatmap(as.matrix(all.list), name='mat', col=colorRamp2(c(0, 1), c("white", "red")),
row_names_side = "left", cluster_rows = FALSE, show_column_dend = FALSE, column_names_rot = 45, width = unit(12, "cm"),
column_title = "Downregulated pathways")
dev.off()

pdf('APR/chrX.down.Hallmark.heatmap.pdf', width=12, height=10)
Heatmap(as.matrix(chrX.list), name='mat', col=colorRamp2(c(0, 1), c("white", "red")),
row_names_side = "left", cluster_rows = FALSE, show_column_dend = FALSE, column_names_rot = 45, width = unit(10, "cm"),
column_title = "Downregulated chrX pathways")
dev.off()






# Identify cell type clusters in chrX genes
cluster <- hclust(dist(plot.matrix[,'pDC']))
cutree(cluster, k=2)

up.chrX <- unlist(lapply(deg, function(x) subset(x, gene %in% rownames(chrX) & logFC.disease_vs_control > 0.5)$gene))

# library(factoextra)
# fviz_nbclust(d, FUNcluster = function(x) cutree(hc, k = x), method = c("wss", "silhouette", "gap_stat"))


x <- strsplit(pathway.analysis$`Tem/Trm cytotoxic T cells`[[1]][2,'geneID'], '/')[[1]]
subset(deg$`Tem/Trm cytotoxic T cells`, gene %in% rownames(chrX)[rownames(chrX) %in% x])
