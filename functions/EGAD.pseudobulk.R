library(EGAD)
library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)

load('/directflow/SCCGGroupShare/projects/lacgra/CoExpNets/bin/run_GBA.Rdata')
load('/directflow/SCCGGroupShare/projects/lacgra/CoExpNets/bin/run_GBA.Rdata')
source('/directflow/SCCGGroupShare/projects/lacgra//CoExpNets/bin/helper_functions.r')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')

# Obtain Phenocarta
phenocarta <- get_phenocarta(species="human", type="all")
biogrid <- get_biogrid(species="9606")

# Build annotations
data("example_annotations")
annotations <- make_annotations(example_annotations$gene2Annot, 
                                example_annotations$genes, 
                                example_annotations$annotationlist)


# Load data and psuedobulk
pbmc <- readRDS('pbmc.female.control-managed.RDS')

cell = 'Tem/Temra cytotoxic T cells'
pbmc.cell <- subset(pbmc, cellTypist == cell)

# Keep genes with expression in 5% of cells
keep <- rowSums(pbmc.cell@assays$RNA@counts > 0) > ncol(pbmc.cell) * 0.05
features <- names(keep[keep == T])
pbmc.cell <- subset(pbmc.cell, features=features)

# Psudobulking by summing counts
expr <- AggregateExpression(pbmc.cell, group.by='individual', slot='counts')$RNA
# expr <- expr[,(colSums(expr) > 0)]
# expr <- expr[rownames(expr) %in% rownames(chrX),]

# Convert gene names to entrez with biomaRt
library(biomaRt)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c('entrezgene_id', 'external_gene_name'), mart=mart)
genes <- genes[!duplicated(genes$external_gene_name),]
genes <- genes[genes$external_gene_name %in% rownames(expr),]
rownames(expr) <- genes[match(rownames(expr), genes$external_gene_name),1]
expr <- expr[!is.na(rownames(expr)),]

# Calculate network edges using correlations and rank standardize
network = EGAD::build_coexp_network(expr, rownames(expr))

# Neighbor Voting
gba_auc_nv <- data.frame(neighbor_voting(annotations, network, nFold=3, output="AUROC"))
gba_auc_nv_flt <- subset(gba_auc_nv, auc > 0.8)
GO_entrez <- lapply(rownames(gba_auc_nv_flt), function(x) subset(example_annotations$gene2Annot, GO %in% x))
names(GO_entrez) <- rownames(gba_auc_nv_flt)
GO_gene <- lapply(GO_entrez, function(x) genes[genes$entrezgene_id %in% x$entrezID,'external_gene_name'])
GO_gene[GO_gene %in% rownames(chrX)]


# Calculate multifunctionality
multifunc_assessment <- calculate_multifunc(annotations)
auc_mf <- auc_multifunc(annotations, multifunc_assessment[,4])
hist <- plot_distribution(auc_mf, xlab="AUROC", med=FALSE, avg=FALSE)

# Node Degree
nd <- node_degree(network)
pdf('EGAD.node.degree.pdf')
hist <- plot_distribution(nd, xlab="Node degree", med = FALSE)
dev.off()

# Network Assortativity
assort <- assortativity(network)
