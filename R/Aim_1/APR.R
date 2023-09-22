library(Seurat)
library(speckle)
library(ggplot2)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(reshape2)
library(Nebulosa)
library(dplyr)

pbmc <- readRDS('pbmc.female.RDS')
#pbmc <- readRDS('pbmc.female.control-managed.RDS')

# Perform propellor cell type abundance testing
output.asin <- propeller(clusters=pbmc$cellTypist, sample=pbmc$individual, group=pbmc$condition, trend=TRUE, transform='asin')
output.logit <- propeller(clusters=pbmc$cellTypist, sample=pbmc$individual, group=pbmc$condition, trend=TRUE, transform='logit')
colnames(output.asin)[1] <- 'celltype'

write.table(output.asin, 'propeller.asin.txt')

# Plot percentage of each celltype in each individual
cell_order <- c(
  "CRTAM+ gamma-delta T cells", "MAIT cells", "Regulatory T cells", 
  "Tcm/Naive cytotoxic T cells", "Tcm/Naive helper T cells", 
  "Tem/Effector helper T cells", "Tem/Temra cytotoxic T cells", 
  "Tem/Trm cytotoxic T cells", "Naive B cells", "Memory B cells", 
  "Plasma cells", "CD16+ NK cells", "NK cells", "Classical monocytes", 
  "Non-classical monocytes", "DC2", "pDC", "Mast cells"
)
# Count number of cellTypist in each individual and divide by total cells in individual
celltype.perc <- pbmc@meta.data[,c('individual', 'condition', 'cellTypist')] %>%
  group_by(individual, cellTypist) %>%
  summarise(n=n()) %>%
  mutate(perc=n/sum(n)*100) %>%
  ungroup() %>%
  data.frame()

# Match celltype.perc$individual to condition
metadata <- unique(pbmc@meta.data[,c('individual', 'condition')])
celltype.perc$condition <- metadata[match(celltype.perc$individual, metadata$individual),'condition']

pdf('APR/cellperc.boxplot.pdf', width=10, height=10)
ggplot(celltype.perc, aes(x=cellTypist, y=perc, fill=condition)) + 
    geom_boxplot() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = unit(c(1, 1, 1, 3), "lines")) + # Adjust the right margin
    labs(x='', y='Percentage of cells (%)', fill='Condition', title='SLE Cell type proportions')
dev.off()

# Load in edgeR results
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
deg <- deg.list('differential.expression/edgeR',logfc=0.5)
deg.chrX <- lapply(deg, function(x) subset(x, gene %in% rownames(chrX)))

DEG <- unique(unlist(lapply(deg, function(x) x$gene)))
DEG.up <- unique(unlist(lapply(deg, function(x) subset(x, logFC.disease_vs_control > 0)$gene)))
DEG.down <- unique(unlist(lapply(deg, function(x) subset(x, logFC.disease_vs_control < 0)$gene)))
DEG.chrX <- unique(unlist(lapply(deg.chrX, function(x) x$gene)))

# Plot heatmap of DEGs
pseudobulk.list <- list()
for(cell in levels(pbmc)){
  tmp <- subset(pbmc, cellTypist==cell)
  tmp <- AggregateExpression(tmp, features=DEG, group.by='individual', slot='counts')$RNA
  pseudobulk.list[[cell]] <- tmp
}

pseudobulk.list[[1]][1:3,1:3]

pdf('APR/DEG.heatmap.pdf', width=10, height=10)


# UpSet plot up upregulated genes
up.lst <- lapply(deg, function(x) subset(x, logFC.disease_vs_control > 0)$gene)
up.mtx <- fromList(up.lst)
pdf('APR/DEG.up.upset.pdf', width=10, height=10)
upset(up.mtx, order.by='freq', nsets=length(up.mtx))
dev.off()

# Heatmap of upregulated genes
pdf('APR/DEG.up.heatmap.pdf', width=10, height=10)
Heatmap(t(up.mtx), col=colorRamp2(c(0, 1), c("white", "red")), show_column_names = FALSE, 
row_names_side = "left", row_dend_side = "right")
dev.off()

# UpSet plot downregulated genes
down.lst <- lapply(deg, function(x) subset(x, logFC.disease_vs_control < 0)$gene)
down.mtx <- fromList(down.lst)
pdf('APR/DEG.down.upset.pdf', width=10, height=10)
upset(down.mtx, order.by='freq', nsets=length(down.mtx))
dev.off()

# Heatmap of downregulated genes
pdf('APR/DEG.down.heatmap.pdf', width=10, height=10)
Heatmap(t(down.mtx), col=colorRamp2(c(0, 1), c("white", "red")), show_column_names = FALSE, 
row_names_side = "left", row_dend_side = "right")
dev.off()

# Read in disgene
disgene <- read.delim('../../DisGeNet/SLE.tsv')

lapply(edgeR, function(x) subset(x, gene == 'TLR7'))

# Gene Set Enrichment analysis
library(fgsea)
hallmark <- gmtPathways('../../gene.sets/h.all.v7.5.1.symbols.gmt')
unique(hallmark$term)

x <- edgeR[[1]]
x <- x[order(x$logFC.disease_vs_control), c('gene', 'logFC.disease_vs_control')]
ranks <- as.list(x$logFC.disease_vs_control)
names(ranks) <- x$gene

fgseaRes <- fgsea(pathways = hallmark, 
                  stats = x,
                  minSize=15,
                  maxSize=500)

gsea.list <- lapply(edgeR, function(x){
  tmp <- x[order(x$logFC.disease_vs_control, decreasing = TRUE),'gene']
  GSEA(tmp, pAdjustMethod='fdr', TERM2GENE=hallmark)
})

x <- x[order(x$logFC.disease_vs_control, decreasing=T), 1:2]
y <- GSEA(x, pAdjustMethod='fdr', TERM2GENE=hallmark)
# Plot density of gene sets for control and disease separately





ISG <- c(
    "ISG15", "IFI6", "IFI44L", "IFI44", "RSAD2", "CXCL10", "IFIT2", "IFIT3", "IFIT1", "IFITM3", "OAS1", "OAS3", "OAS2", "OASL", "EPSTI1", 
    "RNASE1", "RNASE2", "IFI27", "XAF1", "LGALS3BP", "SIGLEC1", "USP18", "APOBEC3A", "APOBEC3B", "MX1"
)

TNFa <- subset(hallmark, term=='HALLMARK_TNFA_SIGNALING_VIA_NFKB')$gene
TGFb <- subset(hallmark, term=='HALLMARK_TGF_BETA_SIGNALING')$gene
INFa <- subset(hallmark, term=='HALLMARK_INTERFERON_ALPHA_RESPONSE')$gene
INGy <- subset(hallmark, term=='HALLMARK_INTERFERON_GAMMA_RESPONSE')$gene


CTL <- c(
  'PRF1', 'GZMH', 'GZMB'
)

control <- subset(pbmc, condition=='control')
disease <- subset(pbmc, condition=='disease')

p1 <- plot_density(control, ISG, joint = TRUE)
pdf('APR/ISG.Nebulosa.control.pdf', width=10, height=10)
p1[[length(ISG)+1]]
dev.off()

p2 <- plot_density(disease, ISG, joint = TRUE)
pdf('APR/ISG.Nebulosa.disease.pdf', width=10, height=10)
p2[[length(ISG)+1]]
dev.off()

pdf('APR/CD19.Nebulosa.pdf')
plot_density(pbmc, 'CD19')
dev.off()



pdf('APR/TLR7.vlnplot.pdf')
VlnPlot(pbmc, features = c('TLR7'), pt.size=0.1, ncol=1, group.by='cellTypist', split.by='condition')
dev.off()



