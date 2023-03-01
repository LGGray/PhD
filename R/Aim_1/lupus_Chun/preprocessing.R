library(Seurat)
library(SeuratDisk)
library(biomaRt)
library(magrittr)
library(ddqcR)
library(DoubletFinder)

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/lupus_Chun/')

ancestry = commandArgs(trailingOnly=TRUE)

#Convert(paste0(ancestry, '.h5ad'), dest='h5seurat', overwrite=F)
pbmc <- LoadH5Seurat(paste0(ancestry, '.h5seurat'), 
                       meta.data=F, assays='RNA', misc=F)

# Add metadata
metadata <- read.delim(paste0(ancestry,'.metadata.csv'), row.names=1, sep=',')

# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
#                    filters = "ensembl_gene_id", 
#                    values = rownames(expr), 
#                    mart = ensembl)
expr <- GetAssayData(pbmc, assay='RNA')
pbmc <- CreateSeuratObject(expr)

pbmc@meta.data <- cbind(pbmc@meta.data, metadata)
pbmc$condition <- ifelse(pbmc$disease == 'systemic lupus erythematosus', 'disease', 'control')

# Remove obvious bad quality cells
pbmc <- initialQC(pbmc)

# Return dataframe of filtering statistics
pdf(paste0('ddqc.plot.', ancestry, '.pdf'))
df.qc <- ddqc.metrics(pbmc)
dev.off()

# Filter out the low quality cells
pbmc <- filterData(pbmc, df.qc)

# Split object by individual for doublet detection
pbmc.split <- SplitObject(pbmc, split.by = "sample_uuid")
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
              project = "lupus_Chun")

# SCTransform
pbmc <- SCTransform(pbmc, verbose = FALSE)

# Cell type clustering and inspection of markers
pbmc <- RunPCA(pbmc)
min.pc <- calc.min.pc(pbmc)
pbmc <- FindNeighbors(pbmc, dims=1:min.pc)
pbmc <- FindClusters(pbmc, resolution=0.5)
pbmc <- RunUMAP(pbmc, dims = 1:min.pc)

pdf(paste0('seurat.clusters.', ancestry, '.DimPlot.pdf'))
DimPlot(pbmc, reduction='umap')
dev.off()

# pbmc.markers <- FindAllMarkers(pbmc, only.pos=T, min.pct=0.25, logfc.threshold = 0.25)
# write.table(pbmc.markers, paste0('FindAllMarkers.', ancestry, '.txt'), row.names=T, quote=F, sep='\t')

mtx <- as.matrix(GetAssayData(pbmc))
write.csv(mtx, paste0('raw.counts.', ancestry, '.csv'))

# Save RDS file for downstream cellTypist analysis
saveRDS(pbmc, paste0('pbmc.', ancestry, '.unlabelled.RDS'))