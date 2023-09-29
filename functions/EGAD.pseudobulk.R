library(EGAD)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

load('/directflow/SCCGGroupShare/projects/lacgra/CoExpNets/bin/run_GBA.Rdata')
source('/directflow/SCCGGroupShare/projects/lacgra/CoExpNets/bin/helper_functions.r')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')

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
assay <- commandArgs(trailingOnly = TRUE)[2]

# Subset for cell type, remove lowly expressed genes, psuedobulk and perform EGAD
result.list <- list()
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
  control.expr <- AggregateExpression(control, group.by='individual', slot='data')[[assay]]
  # Remove genes with stdev = 0
  control.expr <- control.expr[apply(control.expr, 1, sd) > 0,]

  # Disease data
  disease <- subset(pbmc.subset, condition == 'disease')
  # Psudobulking by summing counts
  disease.expr <- AggregateExpression(disease, group.by='individual', slot='data')[[assay]]
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
  result.list[[cell]] <- gba_auc_nv
  # Merge node degree results
  nd.df <- merge(control.nd, disease.nd, by='row.names')
  colnames(nd.df) <- c('gene', 'nd.control', 'nd.disease')
  # Save results
  write.table(gba_auc_nv, paste0('EGAD/', gsub("/|-| ", "_", cell), '.gba.txt'), sep='\t', quote=F)
  write.table(nd.df, paste0('EGAD/', gsub("/|-| ", "_", cell), '.nd.txt'), sep='\t', quote=F)
}

egad.files <- list.files('EGAD', pattern='gba.txt', full.names=T)
egad.result <- lapply(egad.files, read.delim)
names(egad.result) <- gsub('.gba.txt', '', basename(egad.files))

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

disease.pathway <- dplyr::bind_rows(lapply(egad.result, function(x){
  tmp <- subset(x, auc.disease > 0.7 & auc.control < 0.6)
  tmp <- cbind(pathway=rownames(tmp), tmp)
  rownames(tmp) <- NULL
  return(tmp)
  }), .id='celltype')
disease.pathway

chrX.pathway <- lapply(unique(disease.pathway$pathway), function(x){
  features <- subset(hallmark, term == x)$gene
  features[features %in% rownames(chrX)]
})
names(chrX.pathway) <- unique(disease.pathway$pathway)

for(line in 1:nrow(disease.pathway)){
  print(disease.pathway[line,])
  features <- subset(hallmark, term == disease.pathway[line,'pathway'])$gene
  chrX.features <- features[features %in% rownames(chrX)]
  control.coexp <- readRDS(paste0('EGAD/',disease.pathway[line,'celltype'],'.control.network.RDS'))
  control.coexp_gene.set <- control.coexp[colnames(control.coexp) %in% features, rownames(control.coexp) %in% features]
  disease.coexp <- readRDS(paste0('EGAD/',disease.pathway[line,'celltype'],'.disease.network.RDS'))
  disease.coexp_gene.set <- disease.coexp[colnames(disease.coexp) %in% features, rownames(disease.coexp) %in% features]
  if(length(chrX.features) > 0){
    print(chrX.features[chrX.features %in% rownames(control.coexp_gene.set)])
  }
  pdf(paste0('EGAD/',disease.pathway[line,'celltype'],'.', disease.pathway[line,'pathway'], '.control.heatmap.pdf'))
  print(Heatmap(control.coexp_gene.set, name='Spearmans Rho', column_title = "Control", 
  clustering_distance_rows = "spearman", clustering_distance_columns = 'euclidean', 
  clustering_method_rows = "complete", clustering_method_columns = "complete", show_row_names = FALSE, show_column_names = FALSE))
  dev.off()
  pdf(paste0('EGAD/',disease.pathway[line,'celltype'],'.', disease.pathway[line,'pathway'], '.disease.heatmap.pdf'))
  print(Heatmap(disease.coexp_gene.set, name='Spearmans Rho', column_title = "Disease", 
  clustering_distance_rows = "spearman", clustering_distance_columns = 'euclidean',
  clustering_method_rows = "complete", clustering_method_columns = "complete", show_row_names = FALSE, show_column_names = FALSE))
  dev.off()
}