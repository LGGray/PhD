library(EGAD)
library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)
library(biomaRt)

load('/directflow/SCCGGroupShare/projects/lacgra/CoExpNets/bin/run_GBA.Rdata')
load('/directflow/SCCGGroupShare/projects/lacgra/CoExpNets/bin/run_GBA.Rdata')
source('/directflow/SCCGGroupShare/projects/lacgra//CoExpNets/bin/helper_functions.r')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')

# Build annotations
data("example_annotations")
hallmark <- clusterProfiler::read.gmt('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/h.all.v7.5.1.entrez.gmt')
hallmark <- hallmark[,c(2,1)]
annotations <- make_annotations(hallmark, unique(hallmark$gene), unique(hallmark$term))

# annotations <- make_annotations(example_annotations$gene2Annot, 
#                                 example_annotations$genes, 
#                                 example_annotations$annotationlist)

# Calculate multifunctionality
multifunc_assessment <- calculate_multifunc(annotations)
auc_mf <- auc_multifunc(annotations, multifunc_assessment[,4])
hist <- plot_distribution(auc_mf, xlab="AUROC", med=FALSE, avg=FALSE)

# Load data, subset for celltype and condition then psuedobulk
pbmc <- readRDS('pbmc.female.control-managed.RDS')
cell = 'Plasma cells'

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

# Convert gene names to entrez with biomaRt
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c('entrezgene_id', 'external_gene_name'), mart=mart)
genes <- genes[!duplicated(genes$external_gene_name),]
# Control
rownames(control.expr) <- genes[match(rownames(control.expr), genes$external_gene_name),1]
control.expr <- control.expr[!is.na(rownames(control.expr)),]
# Disease
rownames(disease.expr) <- genes[match(rownames(disease.expr), genes$external_gene_name),1]
disease.expr <- disease.expr[!is.na(rownames(disease.expr)),]

### EGAD on Control ###
# Calculate network edges using correlations and rank standardize
control.network = EGAD::build_coexp_network(control.expr, rownames(control.expr))
# Neighbor Voting
control.gba_auc_nv <- data.frame(neighbor_voting(annotations, control.network, nFold=3, output="AUROC"))
# Convert GO ID to term
goterms <- Term(GOTERM)
control.gba_auc_nv$GO_term <- goterms[rownames(control.gba_auc_nv)]
control.gba_auc_nv_flt <- subset(control.gba_auc_nv, auc > 0.6)
control.GO_entrez <- lapply(rownames(control.gba_auc_nv_flt), function(x) subset(example_annotations$gene2Annot, GO %in% x))
names(control.GO_entrez) <- rownames(control.gba_auc_nv_flt)
control.GO_gene <- lapply(control.GO_entrez, function(x) genes[genes$entrezgene_id %in% x$entrezID,'external_gene_name'])

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
# Convert GO ID to term
goterms <- Term(GOTERM)
disease.gba_auc_nv$GO_term <- goterms[rownames(gba_auc_nv)]
disease.gba_auc_nv_flt <- subset(disease.gba_auc_nv, auc > 0.6)
disease.GO_entrez <- lapply(rownames(disease.gba_auc_nv_flt), function(x) subset(example_annotations$gene2Annot, GO %in% x))
names(disease.GO_entrez) <- rownames(disease.gba_auc_nv_flt)
disease.GO_gene <- lapply(disease.GO_entrez, function(x) genes[genes$entrezgene_id %in% x$entrezID,'external_gene_name'])

# Node Degree
disease.nd <- node_degree(disease.network)
pdf('EGAD.node.degree.pdf')
hist <- plot_distribution(nd, xlab="Node degree", med = FALSE)
dev.off()

# Network Assortativity
disease.assort <- assortativity(disease.network)

# setdiff between control and disease
setdiff(rownames(control.gba_auc_nv_flt), rownames(disease.gba_auc_nv_flt))
# intersect between control and disease
intersect(rownames(control.gba_auc_nv_flt), rownames(disease.gba_auc_nv_flt))
