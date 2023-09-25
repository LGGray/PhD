library(EGAD)
library(Seurat)
library(dplyr)
library(ggplot2)

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

# Subset for cell type, remove lowly expressed genes, psuedobulk and perform EGAD
result.list <- list()
for(cell in levels(pbmc)){
  pbmc.subset <- subset(pbmc, cellTypist == cell)
  # Keep genes with expression in 5% of cells
  keep <- rowSums(pbmc.subset@assays$decontXcounts@data > 0) > ncol(pbmc.subset) * 0.05
  features <- names(keep[keep == T])
  pbmc.subset <- subset(pbmc.subset, features=features)

  # Control data
  control <- subset(pbmc.subset, condition == 'control')
  # Psudobulking by summing counts
  control.expr <- AggregateExpression(control, group.by='individual', slot='data')$decontXcounts

  # Disease data
  disease <- subset(pbmc.subset, condition == 'disease')
  # Psudobulking by summing counts
  disease.expr <- AggregateExpression(disease, group.by='individual', slot='data')$decontXcounts

  ### EGAD on Control ###
  # Calculate network edges using correlations and rank standardize
  control.network = EGAD::build_coexp_network(control.expr, rownames(control.expr))
  # Neighbor Voting
  control.gba_auc_nv <- data.frame(neighbor_voting(annotations, control.network, nFold=3, output="AUROC"))
  # Node Degree
  control.nd <- data.frame(node_degree(control.network))

  ## EGAD on Disease ###
  # Calculate network edges using correlations and rank standardize
  disease.network = EGAD::build_coexp_network(disease.expr, gene.list=rownames(disease.expr))
  # Neighbor Voting
  disease.gba_auc_nv <- data.frame(neighbor_voting(annotations, disease.network, nFold=3, output="AUROC"))
  # Node Degree
  disease.nd <- data.frame(node_degree(disease.network))

  # Merge gba results
  gba_auc_nv <- merge(control.gba_auc_nv, disease.gba_auc_nv, by='row.names', suffixes=c('.control', '.disease'))
  rownames(gba_auc_nv) <- gba_auc_nv$Row.names
  gba_auc_nv <- gba_auc_nv[,-1]
  result.list[[cell]] <- gba_auc_nv
  # Merge node degree results
  nd.df <- merge(control.nd, disease.nd, by='row.names') %>% 
    rename(gene=Row.names)
  # Save results
  write.table(gba_auc_nv, paste0('EGAD/', gsub("/|-| ", "_", cell), 'gba.txt'), sep='\t', quote=F)
  write.table(nd.df, paste0('EGAD/', gsub("/|-| ", "_", cell), '.nd.txt'), sep='\t', quote=F)
  # pdf(paste0('EGAD/',gsub("/|-| ", "_", cell), '.gba_auc.scatter.pdf'))
  # ggplot(gba_auc_nv, aes(x=auc.control, y=auc.disease)) + 
  #   geom_point(color=ifelse(gba_auc_nv$auc.disease > 0.8, "red", "black")) +
  #   geom_text(aes(label=ifelse(auc.disease > 0.8, rownames(gba_auc_nv), "")), vjust=-1) +
  #   labs(x="Control AUC", y="Disease AUC") + ggtitle(paste("Guilt by association AUC:", cell))
  # dev.off()
}

# cell <- names(result.list)[15]
# pdf(paste0('EGAD/',gsub("/|-| ", "_", cell), '.gba_auc.scatter.pdf'))
#   ggplot(result.list[[cell]], aes(x=auc.control, y=auc.disease)) + 
#     geom_point(color=ifelse(gba_auc_nv$auc.disease > 0.8, "red", "black")) +
#     geom_text(aes(label=ifelse(auc.disease > 0.8, rownames(gba_auc_nv), "")), vjust=-1, hjust=1) +
#     labs(x="Control AUC", y="Disease AUC") + ggtitle(paste("Guilt by association AUC:", cell))
# dev.off()

# egad.files <- list.files('EGAD', pattern='txt', full.names=T)
# egad.result <- lapply(egad.files, read.delim)
# names(egad.result) <- gsub('.txt', '', basename(egad.files))

# disease.pathway <- dplyr::bind_rows(lapply(egad.result, function(x){
#   tmp <- subset(x, auc.disease > 0.8 & auc.control < 0.8)
#   tmp <- cbind(pathway=rownames(tmp), tmp)
#   rownames(tmp) <- NULL
#   return(tmp)
#   }), .id='celltype')
# disease.pathway

# # identify terms in disease with high AUC
# terms.flt <- rownames(subset(disease.gba_auc_nv, auc > 0.8))
# if(length(terms.flt) > 0){
#   lapply(terms.flt, function(x) {
#     features <- subset(hallmark, term == x)$gene
#     features <- features[features %in% rownames(chrX)]
#     subset(multifunc_assessment, Gene %in% features)
#   })
# }