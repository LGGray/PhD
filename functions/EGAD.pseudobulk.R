library(EGAD)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

load('/directflow/SCCGGroupShare/projects/lacgra/CoExpNets/bin/run_GBA.Rdata')
source('/directflow/SCCGGroupShare/projects/lacgra/CoExpNets/bin/helper_functions.r')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')

# Build annotations
hallmark <- clusterProfiler::read.gmt('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/h.all.v7.5.1.symbols.gmt')
hallmark <- hallmark[,c(2,1)]
annotations <- make_annotations(hallmark, unique(hallmark$gene), unique(hallmark$term))

# # Calculate multifunctionality
# multifunc_assessment <- calculate_multifunc(annotations)
# auc_mf <- auc_multifunc(annotations, multifunc_assessment[,4])
# pdf('multifunc.hist.pdf')
# plot_distribution(auc_mf, xlab="AUROC", med=FALSE, avg=FALSE)
# dev.off()

# Load data, subset for celltype and condition then psuedobulk
pbmc <- readRDS(commandArgs(trailingOnly = TRUE)[1])
assay <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

# Subset for cell type, remove lowly expressed genes, psuedobulk and perform EGAD
for(cell in levels(pbmc)){
  print(cell)
  pbmc.subset <- subset(pbmc, cellTypist == cell)

  if(length(unique(pbmc.subset$condition)) != 2){
    print("Not enough conditions")
    next
  }

  # Keep genes with expression in 5% of cells
  keep <- rowSums(pbmc.subset@assays[[assay]]@data > 0) > ncol(pbmc.subset) * 0.05
  features <- names(keep[keep == T])
  pbmc.subset <- subset(pbmc.subset, features=features)

  # Control data
  control <- subset(pbmc.subset, condition == 'control')
  # Psudobulking by summing counts
  control.expr <- AggregateExpression(control, group.by='individual', slot='counts')[[assay]]
  # Remove genes with stdev = 0
  control.expr <- control.expr[apply(control.expr, 1, sd) > 0,]

  # Disease data
  disease <- subset(pbmc.subset, condition == 'disease')
  # Psudobulking by summing counts
  disease.expr <- AggregateExpression(disease, group.by='individual', slot='counts')[[assay]]
  # Remove genes with stdev = 0
  disease.expr <- disease.expr[apply(disease.expr, 1, sd) > 0,]

  ### EGAD on Control ###
  # Calculate network edges using correlations and rank standardize
  control.network = EGAD::build_coexp_network(control.expr, rownames(control.expr))
  saveRDS(control.network, paste0('EGAD/', gsub("/|-| ", "_", cell), '.control.network.RDS'))
  # Neighbor Voting
  control.gba_auc_nv <- data.frame(neighbor_voting(annotations, control.network, nFold=10, output="AUROC"))
  # Node Degree
  control.nd <- data.frame(node_degree(control.network))

  ## EGAD on Disease ###
  # Calculate network edges using correlations and rank standardize

  disease.network = EGAD::build_coexp_network(disease.expr, gene.list=rownames(disease.expr))
  saveRDS(disease.network, paste0('EGAD/', gsub("/|-| ", "_", cell), '.disease.network.RDS'))
  # Neighbor Voting
  disease.gba_auc_nv <- data.frame(neighbor_voting(annotations, disease.network, nFold=10, output="AUROC"))
  # Node Degree
  disease.nd <- data.frame(node_degree(disease.network))

  # Merge gba results
  gba_auc_nv <- merge(control.gba_auc_nv, disease.gba_auc_nv, by='row.names', suffixes=c('.control', '.disease'))
  rownames(gba_auc_nv) <- gba_auc_nv$Row.names
  gba_auc_nv <- gba_auc_nv[,-1]
  # Merge node degree results
  nd.df <- merge(control.nd, disease.nd, by='row.names')
  colnames(nd.df) <- c('gene', 'nd.control', 'nd.disease')
  # Save results
  write.table(gba_auc_nv, paste0('EGAD/', gsub("/|-| ", "_", cell), '.gba.txt'), sep='\t', quote=F)
  write.table(nd.df, paste0('EGAD/', gsub("/|-| ", "_", cell), '.nd.txt'), sep='\t', quote=F)
}

# Read in result
egad.files <- list.files('EGAD', pattern='gba.txt', full.names=T)
egad.result <- lapply(egad.files, read.delim)
names(egad.result) <- gsub('.gba.txt', '', basename(egad.files))

# plot scatter plot
for(cell in names(egad.result)){
  pdf(paste0('EGAD/', cell, '.scatter.pdf'))
  file <- egad.result[[cell]]
  print(ggplot(file, aes(x=auc.control, y=auc.disease)) + 
  geom_point(color=ifelse(file$auc.disease > 0.7 & file$auc.control < 0.6, "red", "black")) +
  geom_text(aes(label=ifelse(auc.disease > 0.7 & auc.control < 0.6, rownames(file), "")), vjust = "inward", hjust = "inward") +
  geom_hline(yintercept=0.7, linetype="dashed", color="blue") +
  geom_vline(xintercept=0.6, linetype="dashed", color="blue") +
  labs(x="Control AUC", y="Disease AUC") + ggtitle(paste("Guilt by association AUC:", cell)))
  dev.off()
}

# Combine results
disease.pathway <- dplyr::bind_rows(lapply(egad.result, function(x){
  tmp <- subset(x, auc.disease > 0.7 & auc.control < 0.6)
  tmp <- cbind(pathway=rownames(tmp), tmp)
  rownames(tmp) <- NULL
  return(tmp)
  }), .id='celltype')

# print heatmap for control and disease for each pathway
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
for(line in 1:nrow(disease.pathway)){
  features <- subset(hallmark, term == disease.pathway[line,'pathway'])$gene
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

  nd <- read.delim(paste0('EGAD/', cell, '.nd.txt'))
  nd$chrX <- ifelse(nd$gene %in% rownames(chrX), TRUE, FALSE)
  nd$escape <- ifelse(nd$gene %in% rownames(escape), TRUE, FALSE)
  
  features <- subset(hallmark, term == pathway)$gene
  nd.pathway <- subset(nd, gene %in% features)

  # Local Pathway Analysis
  if(sum(nd.pathway$gene %in% rownames(chrX)) == 0){
    print('No chrX genes in pathway')
    local.chrX.test <- NULL
    local.escape.test <- NULL
  } else{
  local.model <- lm(nd.disease ~ nd.control, data = nd.pathway)
  local.residuals <- resid(local.model)
  names(local.residuals) <- nd.pathway$gene 
  local.sorted_residuals <- sort(local.residuals, decreasing = TRUE)  
  local.ranks <- rank(local.sorted_residuals)
  local.gene_info <- data.frame(gene = names(local.sorted_residuals),
  residual = local.sorted_residuals, rank = local.ranks)
  local.gene_info$chrX <- ifelse(local.gene_info$gene %in% rownames(chrX), TRUE, FALSE)
  local.gene_info$escape <- ifelse(local.gene_info$gene %in% rownames(escape), TRUE, FALSE)
  # Mann-Whitney U test to see if chrX genes have higher residuals
  local.chrX.test <- wilcox.test(rank ~ chrX, data = local.gene_info, alternative = "greater", exact=FALSE)
  # Mann-Whitney U test to see if escape genes have higher residuals
  local.escape.test <- wilcox.test(rank ~ escape, data = local.gene_info, alternative = "greater", exact=FALSE)
  }

#   if(!is.null(local.chrX.test)){
#   if(local.chrX.test$p.value < 0.05){
#     pdf('EGAD/node.degree.', cell, '.', pathway, '.pdf')
#     ggplot(nd.pathway, aes(x=nd.control, y=nd.disease, colour=chrX)) +
#       geom_point() +
#       geom_abline(intercept = 0, slope = 1) +
#       xlab('Control Node Degree') + ylab('Disease Node Degree') + ggtitle(paste(cell, pathway, sep=': '))
#     dev.off()
#   }
# }

  # Global Pathway Analysis
  if(sum(nd$gene %in% rownames(chrX)) == 0){
    print('No chrX genes in pathway')
    global.chrX.test <- NULL
    global.escape.test <- NULL
  } else{
  global.model <- lm(nd.disease ~ nd.control, data = nd)
  global.residuals <- resid(global.model)
  names(global.residuals) <- nd$gene
  global.sorted_residuals <- sort(global.residuals, decreasing = TRUE)
  global.ranks <- rank(global.sorted_residuals)
  global.gene_info <- data.frame(gene = names(global.sorted_residuals),
  residual = global.sorted_residuals, rank = global.ranks)
  global.gene_info$chrX <- ifelse(global.gene_info$gene %in% rownames(chrX), TRUE, FALSE)
  global.gene_info$escape <- ifelse(global.gene_info$gene %in% rownames(escape), TRUE, FALSE)
  # Mann-Whitney U test to see if chrX genes have higher residuals
  global.chrX.test <- wilcox.test(rank ~ chrX, data = global.gene_info, alternative = "greater", exact=FALSE)
  # Mann-Whitney U test to see if escape genes have higher residuals
  global.escape.test <- wilcox.test(rank ~ escape, data = global.gene_info, alternative = "greater", exact=FALSE)
  }
  

# if(!is.null(global.chrX.test)){
#   if(global.chrX.test$p.value < 0.05){
#     pdf('EGAD/node.degree.', cell, '.pdf')
#     ggplot(nd, aes(x=nd.control, y=nd.disease, colour=chrX)) +
#       geom_point() +
#       geom_abline(intercept = 0, slope = 1) +
#       xlab('Control Node Degree') + ylab('Disease Node Degree') + ggtitle(paste(cell, pathway, sep=': '))
#     dev.off()
#   }
# }

  # Save results
  if(is.null(local.chrX.test)){
    local.chrX.test <- data.frame(p.value=NA)
  }
  if(is.null(global.chrX.test)){
    global.escape.test <- data.frame(p.value=NA)
  }
  df <- data.frame(celltype=cell, pathway=pathway, 
  local.chrX=local.chrX.test$p.value, global.chrX=global.chrX.test$p.value)
  result_list[[line]] <- df
  # lst <- list(local.chrX.test, local.escape.test, global.chrX.test, global.escape.test)
  # names(lst) <- c('local.chrX.test', 'local.escape.test', 'global.chrX.test', 'global.escape.test')
  # result_list[[line]] <- lst
}
names(result_list) <- disease.pathway$celltype

results <- dplyr::bind_rows(result_list)

write.table(results, 'EGAD/node.degree.chrX.enrichment.txt', sep='\t', quote=F)

cell <- subset(results, global.chrX < 0.05)$celltype
nd <- read.delim(paste0('EGAD/', cell, '.nd.txt'))
nd$chrX <- ifelse(nd$gene %in% rownames(chrX), TRUE, FALSE)
nd.subset <- subset(nd, chrX == TRUE)
pdf(paste0('EGAD/node.degree.', cell, '.chrX.pdf'))
ggplot(nd.subset, aes(x=nd.control, y=nd.disease)) +
geom_point() +
geom_abline(intercept = 0, slope = 1) +
xlab('Control Node Degree') + ylab('Disease Node Degree') + ggtitle(cell)
dev.off()


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

