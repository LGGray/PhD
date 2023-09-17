library(EGAD)
library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)

load('/directflow/SCCGGroupShare/projects/lacgra/CoExpNets/bin/run_GBA.Rdata')
source('/directflow/SCCGGroupShare/projects/lacgra/CoExpNets/bin/helper_functions.r')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')

# Build annotations
data("example_annotations")
hallmark <- clusterProfiler::read.gmt('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/h.all.v7.5.1.symbols.gmt')
hallmark <- hallmark[,c(2,1)]
annotations <- make_annotations(hallmark, unique(hallmark$gene), unique(hallmark$term))

# annotations <- make_annotations(example_annotations$gene2Annot, 
#                                 example_annotations$genes, 
#                                 example_annotations$annotationlist)

# Calculate multifunctionality
multifunc_assessment <- calculate_multifunc(annotations)
auc_mf <- auc_multifunc(annotations, multifunc_assessment[,4])
pdf('multifunc.hist.pdf')
plot_distribution(auc_mf, xlab="AUROC", med=FALSE, avg=FALSE)
dev.off()

# Load data, subset for celltype and condition then psuedobulk
pbmc <- readRDS('pbmc.female.RDS')

control.list <- list()
disease.list <- list()
for(cell in levels(pbmc)){
    # Control data
    control <- subset(pbmc, cellTypist == cell & condition == 'control')
    # Keep genes with expression in 5% of cells
    keep <- rowSums(control@assays$RNA@counts > 0) > ncol(control) * 0.05
    features <- names(keep[keep == T])
    control <- subset(control, features=features)
    # Psudobulking by summing counts
    control.expr <- AggregateExpression(control, group.by='individual', slot='counts')$RNA

    # Disease data
    disease <- subset(pbmc, cellTypist == cell & condition == 'disease')
    # Keep genes with expression in 5% of cells
    keep <- rowSums(disease@assays$RNA@counts > 0) > ncol(disease) * 0.05
    features <- names(keep[keep == T])
    disease <- subset(disease, features=features)
    # Psudobulking by summing counts
    disease.expr <- AggregateExpression(disease, group.by='individual', slot='counts')$RNA

    # Calculate network edges using correlations and rank standardize
    control.network = EGAD::build_coexp_network(control.expr, rownames(control.expr))
    # Neighbor Voting
    control.gba_auc_nv <- data.frame(neighbor_voting(annotations, control.network, nFold=3, output="AUROC"))
    # Calculate network edges using correlations and rank standardize
    disease.network = EGAD::build_coexp_network(disease.expr, gene.list=rownames(disease.expr))
    # Neighbor Voting
    disease.gba_auc_nv <- data.frame(neighbor_voting(annotations, disease.network, nFold=3, output="AUROC"))
}

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
deg <- deg.list('differential.expression/edgeR', logfc=0.5)
names(deg) <- c(
  "CD16+ NK cells", "Classical monocytes", "DC1", "DC2", "MAIT cells", "Mast cells", "Memory B cells", 
  "Naive B cells", "NK cells", "Non-classical monocytes", "pDC", "Plasma cells", "Plasmablasts", "Regulatory T cells", 
  "Tcm/Naive cytotoxic T cells", "Tcm/Naive helper T cells", "Tem/Effector helper T cells", "Tem/Temra cytotoxic T cells", 
  "Tem/Trm cytotoxic T cells"
)

cell = 'Memory B cells'

# Control data
control <- subset(pbmc, cellTypist == cell & condition == 'control', )
# Keep genes with expression in 5% of cells
keep <- rowSums(control@assays$RNA@counts > 0) > ncol(control) * 0.05
features <- names(keep[keep == T])
control <- subset(control, features=features)
# Psudobulking by summing counts
control.expr <- AggregateExpression(control, group.by='individual', slot='counts')$decontXcounts

# Disease data
disease <- subset(pbmc, cellTypist == cell & condition == 'disease', features=deg[[cell]]$gene)
# Keep genes with expression in 5% of cells
keep <- rowSums(disease@assays$RNA@counts > 0) > ncol(disease) * 0.05
features <- names(keep[keep == T])
disease <- subset(disease, features=features)
Psudobulking by summing counts
disease.expr <- AggregateExpression(disease, group.by='individual', slot='counts')$decontXcounts

# # Convert gene names to entrez with biomaRt
# mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# genes <- getBM(attributes=c('entrezgene_id', 'external_gene_name'), mart=mart)
# genes <- genes[!duplicated(genes$external_gene_name),]
# # Control
# rownames(control.expr) <- genes[match(rownames(control.expr), genes$external_gene_name),1]
# control.expr <- control.expr[!is.na(rownames(control.expr)),]
# # Disease
# rownames(disease.expr) <- genes[match(rownames(disease.expr), genes$external_gene_name),1]
# disease.expr <- disease.expr[!is.na(rownames(disease.expr)),]

### EGAD on Control ###
# Calculate network edges using correlations and rank standardize
control.network = EGAD::build_coexp_network(control.expr, rownames(control.expr))
# Neighbor Voting
control.gba_auc_nv <- data.frame(neighbor_voting(annotations, control.network, nFold=3, output="AUROC"))
# Convert GO ID to term
# goterms <- Term(GOTERM)
# control.gba_auc_nv$GO_term <- goterms[rownames(control.gba_auc_nv)]
# control.gba_auc_nv_flt <- subset(control.gba_auc_nv, auc > 0.6)
# control.GO_entrez <- lapply(rownames(control.gba_auc_nv_flt), function(x) subset(example_annotations$gene2Annot, GO %in% x))
# names(control.GO_entrez) <- rownames(control.gba_auc_nv_flt)
# control.GO_gene <- lapply(control.GO_entrez, function(x) genes[genes$entrezgene_id %in% x$entrezID,'external_gene_name'])

# Node Degree
control.nd <- node_degree(control.network)
pdf('EGAD.node.degree.pdf')
hist <- plot_distribution(nd, xlab="Node degree", med = FALSE)
dev.off()

# Network Assortativity
control.assort <- assortativity(control.network)

## EGAD on Disease ###
# Calculate network edges using correlations and rank standardize
disease.network = EGAD::build_coexp_network(disease.expr, gene.list=rownames(disease.expr))

# Neighbor Voting
disease.gba_auc_nv <- data.frame(neighbor_voting(annotations, disease.network, nFold=3, output="AUROC"))
# # Convert GO ID to term
# goterms <- Term(GOTERM)
# disease.gba_auc_nv$GO_term <- goterms[rownames(gba_auc_nv)]
# disease.gba_auc_nv_flt <- subset(disease.gba_auc_nv, auc > 0.6)
# disease.GO_entrez <- lapply(rownames(disease.gba_auc_nv_flt), function(x) subset(example_annotations$gene2Annot, GO %in% x))
# names(disease.GO_entrez) <- rownames(disease.gba_auc_nv_flt)
# disease.GO_gene <- lapply(disease.GO_entrez, function(x) genes[genes$entrezgene_id %in% x$entrezID,'external_gene_name'])

# Node Degree
disease.nd <- node_degree(disease.network)
pdf('EGAD.node.degree.pdf')
hist <- plot_distribution(nd, xlab="Node degree", med = FALSE)
dev.off()

# Network Assortativity
disease.assort <- assortativity(disease.network)

control.nd.chrX <- control.nd[names(control.nd) %in% rownames(chrX)]
disease.nd.chrX <- disease.nd[names(disease.nd) %in% rownames(chrX)]
nd.plot <- data.frame(control.nd.chrX, disease.nd.chrX)

# setdiff between control and disease
setdiff(rownames(control.gba_auc_nv_flt), rownames(disease.gba_auc_nv_flt))
# intersect between control and disease
intersect(rownames(control.gba_auc_nv_flt), rownames(disease.gba_auc_nv_flt))

cor(control.gba_auc_nv$auc, disease.gba_auc_nv$auc, method='spearman')

rownames(subset(control.gba_auc_nv, auc > 0.8))
rownames(subset(disease.gba_auc_nv, auc > 0.8))

# 


# plot heatmap
pdf('APR/test.coexp.heatmap.pdf')
Heatmap(disease.network, name='Spearmans Rho z-score', col=colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
show_row_names = FALSE, show_column_names = FALSE)
dev.off()

library(factoextra)
pdf('APR/nbclust.pdf')
fviz_nbclust(scale(disease.expr), kmeans, method = "gap_stat")
dev.off()

?Heatmap

