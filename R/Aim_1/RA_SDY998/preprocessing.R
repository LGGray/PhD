# Script for building and processing Seurat object from RA-SDY998 data

library(Seurat)
library(SeuratDisk)
library(magrittr)
library(ddqcR)
library(DoubletFinder)
library(parallel)

source('/directflow/SCCGGroupShare/projects/lacgra/R_code/functions/calc.min.pc.R')

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/RA_SDY998')

exprMat <- read.delim('celseq_matrix_ru1_reads.tsv', header=T, row.names = 1)
exprMat[is.na(exprMat)] <- 0
meta <- read.delim('celseq_meta.tsv')
colnames(meta)[3] <- 'individual'
colnames(meta)[5] <- 'condition'
meta$condition <- ifelse(meta$condition == 'RA', 'disease', 'control')
pbmc <- CreateSeuratObject(counts = exprMat)
pbmc@meta.data <- cbind(pbmc@meta.data, meta)

# Remove obvious bad quality cells
pbmc <- initialQC(pbmc)

# Return dataframe of filtering statistics
pdf('ddqc.plot.pdf')
df.qc <- ddqc.metrics(pbmc)
dev.off()

# Filter out the cells
pbmc <- filterData(pbmc, df.qc)

pbmc.split <- SplitObject(pbmc, split.by = "individual")

# loop through samples to run DoubletFinder on each individual
for (i in 1:length(pbmc.split)) {
  # print the sample we are on
  print(paste0("Sample ",i))
  
  # Pre-process seurat object with standard seurat workflow
  pbmc.sample <- NormalizeData(pbmc.split[[i]])
  pbmc.sample <- FindVariableFeatures(pbmc.sample)
  pbmc.sample <- ScaleData(pbmc.sample)
  if(ncol(pbmc.sample) > 50){
    pbmc.sample <- RunPCA(pbmc.sample)
  } else{
    pbmc.split[[i]] <- NULL
    break
  }
  
  # Find significant PCs
  stdv <- pbmc.sample[["pca"]]@stdev
  sum.stdv <- sum(pbmc.sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  print(min.pc)
  
  # finish pre-processing
  pbmc.sample <- RunUMAP(pbmc.sample, dims = 1:min.pc)
  pbmc.sample <- FindNeighbors(object = pbmc.sample, dims = 1:min.pc)              
  pbmc.sample <- FindClusters(object = pbmc.sample, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep_v3(pbmc.sample, PCs = 1:min.pc, num.cores = detectCores()/2)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- pbmc.sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(pbmc.sample@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  pbmc.sample <- doubletFinder_v3(seu = pbmc.sample, 
                                   PCs = 1:min.pc, 
                                   pK = optimal.pk,
                                   nExp = nExp.poi.adj)
  metadata <- pbmc.sample@meta.data
  colnames(metadata)[ncol(metadata)] <- "doublet_finder"
  pbmc.sample@meta.data <- metadata 
  
  # subset and save
  pbmc.singlets <- subset(pbmc.sample, doublet_finder == "Singlet")
  pbmc.split[[i]] <- pbmc.singlets
  remove(pbmc.singlets)
}

pbmc <- merge(x = pbmc.split[[1]],
              y = c(pbmc.split[-1]),
              project = "RA_SDY998")

# SCTransform
pbmc <- SCTransform(pbmc, verbose = FALSE)

# Cell type clustering and inspection of markers
pbmc <- RunPCA(pbmc)
min.pc <- calc.min.pc(pbmc)
pbmc <- FindNeighbors(pbmc, dims=1:min.pc)
pbmc <- FindClusters(pbmc, resolution=0.5)
pbmc <- RunUMAP(pbmc, dims = 1:min.pc)

pdf('seurat.clusters.DimPlot.pdf')
DimPlot(pbmc, reduction='umap')
dev.off()

pbmc.markers <- FindAllMarkers(pbmc, only.pos=T, min.pct=0.25, logfc.threshold = 0.25)
write.table(pbmc.markers, 'FindAllMarkers.txt', row.names=T, quote=F, sep='\t')

# Export .h5ad file for cellTypist
SaveH5Seurat(pbmc, filename = "pbmc.h5Seurat", overwrite = TRUE)
Convert("pbmc.h5Seurat", dest = "h5ad", overwrite = TRUE)

mtx <- as.matrix(GetAssayData(pbmc))
write.csv(mtx, 'raw.counts.csv')

# Save RDS file for downstream cellTypist analysis
saveRDS(pbmc, 'pbmc.unlabelled.RDS')