# Script for building and processing Seurat object from UC_GSE182270 data

# Prevent warnings from printing
options(warn=-1)
# Load Libraries
library(Seurat)
library(SeuratDisk)
library(magrittr)
library(ddqcR)
library(DoubletFinder)
library(parallel)

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrY.Rdata')
source('/directflow/SCCGGroupShare/projects/lacgra/R_code/functions/calc.min.pc.R')

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/UC_GSE182270')

files <- list.files(recursive = T, pattern='.tsv|mtx')
files <- files[grep('Raw', files)]

individuals <- dirname(files) %>% 
  gsub('_Grch38/Raw', '', .) %>%
  unique()

expr.list <- lapply(individuals, function(x){
  tmp <- unique(dirname(files[grep(x, files)]))
  Read10X(tmp)
})

pbmc.list <- list()
for(i in 1:length(expr.list)){
  pbmc.list[i] <- CreateSeuratObject(counts = expr.list[[i]])
}
pbmc <- merge(pbmc.list[[1]], c(pbmc.list[[2]], pbmc.list[[3]], pbmc.list[[4]],
                              pbmc.list[[5]], pbmc.list[[6]], pbmc.list[[7]],
                              pbmc.list[[8]]),
              add.cell.ids = individuals)
pbmc$individual <- gsub('_[ACTG]+-.', '', colnames(pbmc))
pbmc$condition <- gsub('\\d.+', '', pbmc$individual)
pbmc$condition <- ifelse(pbmc$condition == 'UC', 'Disease', 'Control')

# Remove obvious bad quality cells
pbmc <- initialQC(pbmc)

# Return dataframe of filtering statistics
pdf('ddqc.plot.pdf')
df.qc <- ddqc.metrics(pbmc)
dev.off()

# Filter out the low quality cells
pbmc <- filterData(pbmc, df.qc)

# Split object by individual for doublet detection
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
              project = "UC_GSE182270")

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