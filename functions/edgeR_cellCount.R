library(edgeR)
library(Seurat)
library(MAST)
library(qvalue)
library(tidyverse)

if(dir.exists('differential.expression') != TRUE){dir.create('differential.expression')}
if(dir.exists('differential.expression/edgeR') != TRUE){dir.create('differential.expression/edgeR')}
if(dir.exists('differential.expression/MAST') != TRUE){dir.create('differential.expression/MAST')}

# Read in file from command line
pbmc <- readRDS(commandArgs(trailingOnly = TRUE)[1])


for (cell in levels(pbmc)){
  # Select cell type
  print(cell)
  # subset object by cell type
  pbmc.cell <- subset(pbmc, cellTypist == cell)
  # check if there are enough cells and skip if not
  if(nrow(pbmc.cell) < 30){
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
  expr <- AggregateExpression(pbmc.cell, group.by='individual', slot='counts')$RNA
  expr <- expr[,(colSums(expr) > 0)]

  # edgeR-QLFTest
  targets = unique(data.frame(condition = pbmc.cell$condition,
                      individual = pbmc.cell$individual))
  targets$cellCount <- cellCount[,cell]
  targets <- targets[match(colnames(expr), targets$individual),]

  # targets$SV1 <- tapply(pbmc.cell$SV1, pbmc.cell$individual, sum)
  # targets$SV2 <- tapply(pbmc.cell$SV2, pbmc.cell$individual, sum)
  rownames(targets) <- targets$individual
  design <- model.matrix(~0 + cellCount + condition, data=targets)
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
      res$FDR <- qvalue(p = res$PValue)$qvalues
  cell = gsub("/|-| ", "_", cell)
  write.table(res, paste0("differential.expression/edgeR/", cell, ".txt"),
              row.names=F, sep="\t", quote = F)
}

print("Done with edgeR-QLF")

deg <- deg.list('differential.expression/edgeR', logfc=0.05)
deg.chrX <- lapply(deg, function(x) subset(x, gene %in% rownames(chrX)))

# Heatmap of chrX genes across celltypes
genes <- unique(unlist(lapply(deg.chrX, function(x) x$gene)))
plot.matrix <- matrix(0, nrow=length(genes), ncol=length(deg.chrX))
rownames(plot.matrix) <- genes
colnames(plot.matrix) <- names(deg.chrX)
# Match genes to rownames
for (i in 1:length(deg.chrX)){
  plot.matrix[match(deg.chrX[[i]]$gene, genes),i] <- deg.chrX[[i]]$logFC.disease_vs_control
}
# Plot heatmap
library(gplots)
pdf('differential.expression/chrX.heatmap.pdf')
heatmap.2(plot.matrix, trace='none', col=rev(colorRampPalette(c('blue', 'white', 'red'))(100)), 
          scale='row', key=F, cexRow=0.5, srtCol=45, cexCol=0.5,  margins=c(10,10))    
dev.off()

edgeR <- deg.list('differential.expression/edgeR', filter=F)
# replace list colnames[2] with 'logFC'
edgeR <- lapply(edgeR, function(x) {colnames(x)[2] <- 'logFC'; return(x)})

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/chisq.test.degs.R')
lapply(edgeR, function(x) chisq.test.edgeR(x, rownames(chrX), 0.05))

# Identify cell type clusters in chrX genes
cluster <- hclust(dist(plot.matrix[,'pDC']))
cutree(cluster, k=2)

# library(factoextra)
# fviz_nbclust(d, FUNcluster = function(x) cutree(hc, k = x), method = c("wss", "silhouette", "gap_stat"))