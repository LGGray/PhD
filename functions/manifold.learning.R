library(Seurat)

pbmc <- readRDS('pbmc.female.RDS')

Memory <- subset(pbmc, idents = "Memory B cells")
Plasma <- subset(pbmc, idents = "Plasma cells")
ABC <- subset(pbmc, idents = "Age-assocciated B cells")

Memory <- RunTSNE(Memory, dims = 1:20)
Plasma <- RunTSNE(Plasma, dims = 1:20)

Memory <- FindNeighbors(Memory, dims = 1:20)
Memory <- FindClusters(Memory, resolution = 0.5)

Plasma <- FindNeighbors(Plasma, dims = 1:20)
Plasma <- FindClusters(Plasma, resolution = 0.5)

# Find genes that are differentially expressed between clusters
Memory.markers <- FindAllMarkers(Memory, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Memory.markers <- subset(Memory.markers, p_val_adj < 0.05)
Plasma.markers <- FindAllMarkers(Plasma, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Plasma.markers <- subset(Plasma.markers, p_val_adj < 0.05)

# find IGH genes in markers
Memory.markers[grep('^IGH[DAMGE]',Memory.markers$gene),]
Plasma.markers[grep('^IGH[DAMGE]',Plasma.markers$gene),]

# Subset Memory by heavy chain expression and sum expression
heavy.chain <- rownames(Memory)[grep('^IGH[DAMGE]',rownames(Memory))]
Memory.heavy <- rowSums(subset(Memory, features = heavy.chain))
# Calculate expression of heavy chain genes in each cluster
AverageExpression(Memory, features = heavy.chain)