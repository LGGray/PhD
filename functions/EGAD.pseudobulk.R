library(EGAD)
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)

load('/directflow/SCCGGroupShare/projects/lacgra/CoExpNets/bin/run_GBA.Rdata')
source('/directflow/SCCGGroupShare/projects/lacgra/CoExpNets/bin/helper_functions.r')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')


# Build annotations - Hallmark
hallmark <- clusterProfiler::read.gmt('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/h.all.v7.5.1.symbols.gmt')
hallmark <- hallmark[,c(2,1)]
hallmark_annotations <- make_annotations(hallmark, unique(hallmark$gene), unique(hallmark$term))

GOBP <- clusterProfiler::read.gmt('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/c5.go.bp.v7.5.1.symbols.gmt')
GOBP <- GOBP[,c(2,1)]
GOBP_annotations <- make_annotations(GOBP, unique(GOBP$gene), unique(GOBP$term))

# # Calculate multifunctionality
multifunc_assessment <- calculate_multifunc(GOBP_annotations)
colnames(multifunc_assessment)[1] <- 'gene'
auc_mf <- auc_multifunc(GOBP_annotations, multifunc_assessment[,4])
# pdf('multifunc.hist.pdf')
# plot_distribution(auc_mf, xlab="AUROC", med=FALSE, avg=FALSE)
# dev.off()

# Load data, subset for celltype and condition then psuedobulk
pbmc <- readRDS(commandArgs(trailingOnly = TRUE)[1])
assay <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

# Subset for cell type, remove lowly expressed genes, psuedobulk and perform EGAD
# for(cell in levels(pbmc)){
#   print(cell)
#   pbmc.subset <- subset(pbmc, cellTypist == cell)

#   if(length(unique(pbmc.subset$condition)) != 2){
#     print("Not enough conditions")
#     next
#   }

#   # Keep genes with expression in 5% of cells
#   keep <- rowSums(pbmc.subset@assays[[assay]]@data > 0) > ncol(pbmc.subset) * 0.05
#   features <- names(keep[keep == T])
#   pbmc.subset <- subset(pbmc.subset, features=features)

#   # Control data
#   control <- subset(pbmc.subset, condition == 'control')
#   # Psudobulking by summing counts
#   control.expr <- AggregateExpression(control, group.by='individual', slot='counts')[[assay]]
#   # Remove genes with stdev = 0
#   control.expr <- control.expr[apply(control.expr, 1, sd) > 0,]

#   # Disease data
#   disease <- subset(pbmc.subset, condition == 'disease')
#   # Psudobulking by summing counts
#   disease.expr <- AggregateExpression(disease, group.by='individual', slot='counts')[[assay]]
#   # Remove genes with stdev = 0
#   disease.expr <- disease.expr[apply(disease.expr, 1, sd) > 0,]
# }

psuedobulked <- list.files('EGAD', pattern='network.RDS', full.names=T)
celltypes <- unique(gsub('\\..+.network.RDS', '', basename(psuedobulked)))
gene.set <- 'GOBP'
for(cell in celltypes){
  print(cell)
  control <- readRDS(paste0('EGAD/', cell, '.control.network.RDS'))
  disease <- readRDS(paste0('EGAD/', cell, '.disease.network.RDS'))
  # Control
  if(gene.set == 'GOBP'){
    control.gba_auc_nv <- data.frame(neighbor_voting(GOBP_annotations, control, nFold=10, output="AUROC"))
  } else{
    control.gba_auc_nv <- data.frame(neighbor_voting(hallmark_annotations, control, nFold=10, output="AUROC"))
  }
  # Node Degree
  control.nd <- data.frame(node_degree(control))
  control.nd$std.nd <- apply(control.nd, 1, function(x) x/nrow(control.nd))

  # Disease
  if(gene.set == 'GOBP'){
    disease.gba_auc_nv <- data.frame(neighbor_voting(GOBP_annotations, disease, nFold=10, output="AUROC"))
  } else{
    disease.gba_auc_nv <- data.frame(neighbor_voting(hallmark_annotations, disease, nFold=10, output="AUROC"))
  }
  # Node Degree
  disease.nd <- data.frame(node_degree(disease))
  disease.nd$std.nd <- apply(disease.nd, 1, function(x) x/nrow(disease.nd))

  # Merge gba results
  gba_auc_nv <- merge(control.gba_auc_nv, disease.gba_auc_nv, by='row.names', suffixes=c('.control', '.disease'))
  rownames(gba_auc_nv) <- gba_auc_nv$Row.names
  gba_auc_nv <- gba_auc_nv[,-1]
  # Merge node degree results
  nd.df <- merge(control.nd, disease.nd, by='row.names')
  colnames(nd.df) <- c('gene', 'nd.control', 'std.nd.control', 'nd.disease', 'std.nd.disease')
  # Save results
  write.table(gba_auc_nv, paste0('EGAD/', gene.set, '/', gsub("/|-| ", "_", cell), '.gba.txt'), sep='\t', quote=F)
  write.table(nd.df, paste0('EGAD/', gene.set, '/', gsub("/|-| ", "_", cell), '.nd.txt'), sep='\t', quote=F)
}

node_degree <- melt(nd.df[,c(1,3,5)], id.vars='gene')
node_degree$variable <- factor(node_degree$variable, levels=c('std.nd.control', 'std.nd.disease'))
pdf(paste0('EGAD/', gene.set, '/', 'node.degree.density.pdf'))
print(ggplot(node_degree, aes(x=value, fill=variable)) + geom_density(alpha=0.5) + facet_wrap(~variable))
dev.off()

pdf(paste0('EGAD/', gene.set, '/', 'node.degree.violin.pdf'))
ggplot(node_degree, aes(y=value, x=variable)) + geom_violin()
dev.off()
   +
  geom_text(data = subset(node_degree, abs(value - mean(value)) > 3 * sd(value)),
            aes(label = gene), vjust = -0.5, hjust = 0.5, color = "red")


wilcox.test(nd.df$std.nd.disease, nd.df$std.nd.control)

# Read in result
egad.files <- list.files(paste0('EGAD/', gene.set), pattern='gba.txt', full.names=T)
egad.result <- lapply(egad.files, read.delim)
names(egad.result) <- gsub('.gba.txt', '', basename(egad.files))


# plot scatter plot
for(cell in names(egad.result)){
  pdf(paste0('EGAD/', gene.set, '/', cell, '.scatter.pdf'))
  file <- egad.result[[cell]]
  print(ggplot(file, aes(x=auc.control, y=auc.disease)) + 
  geom_point(color=ifelse(file$auc.disease > 0.9 & file$auc.control < 0.5, "red", "black")) +
  # geom_text(aes(label=ifelse(auc.disease > 0.9 & auc.control < 0.5, rownames(file), "")), vjust = "inward", hjust = "inward") +
  geom_hline(yintercept=0.9, linetype="dashed", color="blue") +
  geom_vline(xintercept=0.5, linetype="dashed", color="blue") +
  labs(x="Control AUC", y="Disease AUC") + ggtitle(paste("Guilt by association AUC:", cell)))
  dev.off()
}

# Combine results
disease.pathway <- dplyr::bind_rows(lapply(egad.result, function(x){
  tmp <- subset(x, auc.disease > 0.9 & auc.control < 0.5)
  tmp <- cbind(pathway=rownames(tmp), tmp)
  rownames(tmp) <- NULL
  return(tmp)
  }), .id='celltype')

# print heatmap for control and disease for each pathway
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
for(line in 1:nrow(disease.pathway)){
  features <- subset(gene.set, term == disease.pathway[line,'pathway'])$gene
  # chrX.features <- features[features %in% rownames(chrX)]
  control.coexp <- readRDS(paste0('EGAD/',disease.pathway[line,'celltype'],'.control.network.RDS'))
  control.coexp_gene.set <- control.coexp[colnames(control.coexp) %in% features, rownames(control.coexp) %in% features]
  disease.coexp <- readRDS(paste0('EGAD/',disease.pathway[line,'celltype'],'.disease.network.RDS'))
  disease.coexp_gene.set <- disease.coexp[colnames(disease.coexp) %in% features, rownames(disease.coexp) %in% features]

  pdf(paste0('EGAD/',disease.pathway[line,'celltype'],'.', disease.pathway[line,'pathway'], '.heatmap.pdf'), width=10, height=10)
  control_heatmap <- Heatmap(control.coexp_gene.set, column_title = "Control", col=col_fun,
  clustering_distance_rows = "euclidean", clustering_distance_columns = 'euclidean', 
  clustering_method_rows = "complete", clustering_method_columns = "complete", show_heatmap_legend = FALSE)
  disease_heatmap <- Heatmap(disease.coexp_gene.set, column_title = "Disease", col=col_fun,
  clustering_distance_rows = "euclidean", clustering_distance_columns = 'euclidean',
  clustering_method_rows = "complete", clustering_method_columns = "complete", show_heatmap_legend = FALSE)
  pushViewport(viewport(x = 0, width = 0.5, just = "left"))
  draw(control_heatmap, newpage = FALSE)
  popViewport()
  pushViewport(viewport(x = 1, width = 0.5, just = "right"))
  draw(disease_heatmap, newpage = FALSE)
  popViewport()
  dev.off()
}

# features <- subset(hallmark, term == disease.pathway[line,'pathway'])$gene
# chrX.features <- features[features %in% rownames(chrX)]
# control.coexp <- readRDS(paste0('EGAD/',disease.pathway[line,'celltype'],'.control.network.RDS'))
# control.coexp_gene.set <- control.coexp[colnames(control.coexp) %in% features, rownames(control.coexp) %in% features]
# disease.coexp <- readRDS(paste0('EGAD/',disease.pathway[line,'celltype'],'.disease.network.RDS'))
# disease.coexp_gene.set <- disease.coexp[colnames(disease.coexp) %in% features, rownames(disease.coexp) %in% features]
# pdf(paste0('EGAD/',disease.pathway[line,'celltype'],'.', disease.pathway[line,'pathway'], '.heatmap.pdf'), width=10, height=10)
# control_heatmap <- Heatmap(control.coexp_gene.set, column_title = "Control", col=col_fun,
# clustering_distance_rows = "euclidean", clustering_distance_columns = 'euclidean', 
# clustering_method_rows = "complete", clustering_method_columns = "complete", show_heatmap_legend = FALSE, 
# show_row_names = FALSE, show_column_names = FALSE)
# disease_heatmap <- Heatmap(disease.coexp_gene.set, column_title = "Disease", col=col_fun,
# clustering_distance_rows = "euclidean", clustering_distance_columns = 'euclidean',
# clustering_method_rows = "complete", clustering_method_columns = "complete", show_heatmap_legend = FALSE, 
# show_row_names = FALSE, show_column_names = FALSE)
# pushViewport(viewport(x = 0, width = 0.5, just = "left"))
# draw(control_heatmap, newpage = FALSE)
# popViewport()
# pushViewport(viewport(x = 1, width = 0.5, just = "right"))
# draw(disease_heatmap, newpage = FALSE)
# popViewport()
# dev.off()

# Node degree analysis
result_list <- list()
for(line in 1:nrow(disease.pathway)){
  cell <- disease.pathway[line,'celltype']
  pathway <- disease.pathway[line,'pathway']
  print(paste(cell, pathway, sep=': '))

  nd <- read.delim(paste0('EGAD/', gene.set, '/', cell, '.nd.txt'))
  nd$chrX <- ifelse(nd$gene %in% rownames(chrX), TRUE, FALSE)
  
  features <- subset(get(gene.set), term == pathway)$gene
  nd.pathway <- subset(nd, gene %in% features)

  # Local Pathway Analysis
  if(sum(nd.pathway$gene %in% rownames(chrX)) == 0){
    print('No chrX genes in pathway')
    local.chrX.test <- NULL
  } else{
  local.residuals <- nd.pathway$nd.disease - nd.pathway$nd.control
  names(local.residuals) <- nd.pathway$gene 
  local.sorted_residuals <- sort(local.residuals, decreasing = TRUE)  
  local.ranks <- rank(local.sorted_residuals)
  local.gene_info <- data.frame(gene = names(local.sorted_residuals),
  residual = local.sorted_residuals, rank = local.ranks)
  local.gene_info$chrX <- factor(ifelse(local.gene_info$gene %in% rownames(chrX), TRUE, FALSE))
  # Mann-Whitney U test to see if chrX genes have higher residuals
  local.chrX.test <- wilcox.test(residual ~ chrX, data = local.gene_info, exact=FALSE)
  }

  # Global Pathway Analysis
  if(sum(nd$gene %in% rownames(chrX)) == 0){
    print('No chrX genes in pathway')
    global.chrX.test <- NULL
  } else{
  global.residuals <- nd$nd.disease - nd$nd.control
  names(global.residuals) <- nd$gene
  global.sorted_residuals <- sort(global.residuals, decreasing = TRUE)
  global.ranks <- rank(global.sorted_residuals)
  global.gene_info <- data.frame(gene = names(global.sorted_residuals),
  residual = global.sorted_residuals, ranks = global.ranks)
  global.gene_info$chrX <- factor(ifelse(global.gene_info$gene %in% rownames(chrX), TRUE, FALSE))
  # Mann-Whitney U test to see if chrX genes have higher residuals
  global.chrX.test <- wilcox.test(residual ~ chrX, data = global.gene_info, alternative = "greater", exact=FALSE)
  }

  # Save results
  if(is.null(local.chrX.test)){
    local.chrX.test <- data.frame(p.value=NA)
  }
  df <- data.frame(celltype=cell, pathway=pathway, 
  local.chrX=local.chrX.test$p.value, global.chrX=global.chrX.test$p.value)
  result_list[[line]] <- df
}

results <- dplyr::bind_rows(result_list)
write.table(results, paste0('EGAD/', gene.set, '/', 'node.degree.chrX.enrichment.txt'), sep='\t', quote=F)

# Pathways heatmap
pathway_list <- split(disease.pathway, disease.pathway$celltype)
pathway_list <- lapply(pathway_list, function(x) x$pathway)
binary_matrix <- UpSetR::fromList(pathway_list)
rownames(binary_matrix) <- unique(unlist(pathway_list))
pdf(paste0('EGAD/', gene.set, '/', 'heatmap.pdf'))
col_fun = colorRamp2(c(0, 1), c("white", "red"))
Heatmap(as.matrix(binary_matrix), show_column_names=T, show_row_names = FALSE,
col=col_fun)
dev.off()


for(cell in celltypes){
  nd <- read.delim(paste0('EGAD/', gene.set, '/', cell, '.nd.txt'))
  nd.subset <- subset(nd, gene %in% rownames(chrX))
  nd.subset <- merge(nd.subset, multifunc_assessment, by='gene')
  pdf(paste0('EGAD/', gene.set, '/', 'node.degree.', cell, '.chrX.pdf'))
  top5 <- head(nd.subset[order(nd.subset$nd.disease - nd.subset$nd.control, decreasing = TRUE), ], 5) 
  bottom5 <- head(nd.subset[order(nd.subset$nd.disease - nd.subset$nd.control, decreasing = FALSE), ], 5)
  print(ggplot(nd.subset, aes(x = nd.control, y = nd.disease)) +
    geom_point(aes(color=ifelse(gene %in% c(top5$gene, bottom5$gene), 'red', 'black'))) +
    geom_text_repel(data = top5, aes(label = gene), vjust = "inward", hjust = "inward") +
    geom_text_repel(data = bottom5, aes(label = gene), vjust = "inward", hjust = "inward") +
    geom_abline(intercept = 0, slope = 1) +
    theme(legend.position = "none") +
    xlab('Control Node Degree') + ylab('Disease Node Degree') + ggtitle(cell))
  dev.off()
}

for(cell in celltypes){
  print(cell)
  nd <- read.delim(paste0('EGAD/', gene.set, '/', cell, '.nd.txt'))
  nd.subset <- subset(nd, gene %in% rownames(chrX))
  top5 <- head(nd.subset[order(nd.subset$nd.disease - nd.subset$nd.control, decreasing = TRUE), ], 5) 
  bottom5 <- head(nd.subset[order(nd.subset$nd.disease - nd.subset$nd.control, decreasing = FALSE), ], 5)
  print(rownames(escape)[rownames(escape) %in% c(top5$gene, bottom5$gene)])
}

# Plot heatmap 
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
chrX.match <- rownames(chrX)[rownames(chrX) %in% rownames(network)]
pdf('EGAD/Plasmablast.heatmap.pdf')
ha = rowAnnotation(genes = anno_mark(at = match(chrX.match, rownames(network)), 
    labels = chrX.match, labels_gp = gpar(fontsize = 10)))
Heatmap(network, column_title = "Disease: Plasmablasts", col=col_fun, right_annotation = ha,
clustering_distance_rows = "euclidean", clustering_distance_columns = 'euclidean',
clustering_method_rows = "complete", clustering_method_columns = "complete", show_heatmap_legend = TRUE, 
show_row_names = FALSE, show_column_names = FALSE)
dev.off()

