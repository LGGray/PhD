library(GEOquery)
library(Seurat)
library(SeuratDisk)
library(magrittr)
library(ddqcR)
library(DoubletFinder)
library(parallel)

load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrY.Rdata')
source('/directflow/SCCGGroupShare/projects/lacgra/R_code/functions/calc.min.pc.R')

#options(download.file.method.GEOquery = "wget")

#eList2 <- getGEOSuppFiles('GSE134809')
#untar('GSE134809_RAW.tar')

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/CD_GSE134809')

samples <- '_67|_126|_127|_134|_157|_179|_185|_191|_194'
genes.path <- list.files('.', pattern = 'genes.tsv.gz', recursive = T)
genes.path <- genes.path[grep(samples, genes.path)]
matrix.path <- list.files('.', pattern = 'matrix.mtx.gz', recursive = T)
matrix.path <- matrix.path[grep(samples, matrix.path)]
barcodes.path <- list.files('.', pattern='barcodes.tsv.gz', recursive = T)
barcodes.path <- barcodes.path[grep(samples, barcodes.path)]

exp.mtx <- ReadMtx(mtx = matrix.path[1], 
                      features=genes.path[1], 
                      cells=barcodes.path[1])
GSM4761136 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[2], 
                   features=genes.path[2], 
                   cells=barcodes.path[2])
GSM4761137 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[3], 
                   features=genes.path[3], 
                   cells=barcodes.path[3])
GSM4761138 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[4], 
                   features=genes.path[4], 
                   cells=barcodes.path[4])
GSM4761139 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[5], 
                   features=genes.path[5], 
                   cells=barcodes.path[5])
GSM4761140 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[6], 
                   features=genes.path[6], 
                   cells=barcodes.path[6])
GSM4761141 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[7], 
                   features=genes.path[7], 
                   cells=barcodes.path[7])
GSM4761142 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[8], 
                   features=genes.path[8], 
                   cells=barcodes.path[8])
GSM4761143 <- CreateSeuratObject(counts=exp.mtx)
exp.mtx <- ReadMtx(mtx = matrix.path[9], 
                   features=genes.path[9], 
                   cells=barcodes.path[9])
GSM4761144 <- CreateSeuratObject(counts=exp.mtx)

# Merge files together
pbmc <- merge(GSM4761136, y = c(GSM4761137, GSM4761138, GSM4761139,
                                GSM4761140, GSM4761141, GSM4761142,
                                GSM4761143, GSM4761144), 
              add.cell.ids = c('GSM4761136', 'GSM4761137', 'GSM4761138', 
                               'GSM4761139','GSM4761140', 'GSM4761141', 
                               'GSM4761142','GSM4761143', 'GSM4761144'), project = "CD")
# Add individual  and condition to metadata
pbmc$individual <- gsub('_.+', '', rownames(pbmc@meta.data))
pbmc$condition <- 'CD'

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
  sweep.list <- paramSweep_v3(pbmc.sample, PCs = 1:min.pc, num.cores = detectCores() - 1)
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
              project = "CD_GSE134809")

# Export .h5ad file for cellTypist
SaveH5Seurat(pbmc, filename = "pbmc.h5Seurat")
Convert("pbmc.h5Seurat", dest = "h5ad")

mtx <- as.matrix(GetAssayData(pbmc))
write.csv(mtx, 'raw.counts.csv')

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

# Read in PBMC reference dataset
reference <- LoadH5Seurat("/directflow/SCCGGroupShare/projects/lacgra/azimuth.reference/pbmc_multimodal.h5seurat")
# Find anchors between reference and query
anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)
# Transfer cell type labels and protein data from reference to query
pbmc <- MapQuery(
  anchorset = anchors,
  query = pbmc,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

# scPred cell type annotation
reference <- getFeatureSpace(reference, "celltype.l2")
reference <- trainModel(reference)
get_scpred(reference)
pbmc <- scPredict(pbmc, reference)
crossTab(query, 'predicted.celltype.l2', 'scpred_prediction')

Idents(pbmc) <- 'predicted.celltype.l2'

pdf('DimPlot.all.pdf')
DimPlot(pbmc, label = TRUE, reduction='ref.umap', repel=T) + NoLegend()
dev.off()

# Determine sex of individuals by expression of chrY and XIST
set.seed(42)
exp <- AverageExpression(pbmc, assays='SCT', features=c('XIST', rownames(chrY)), group.by='individual')
pca <- prcomp(t(exp[[1]]), scale=T)
meta <- unique(pbmc@meta.data[,c('individual', 'condition')])
pdf('sex.PCA.pdf')
ggplot(data.frame(pca$x), aes(x=PC1, y=PC2, colour=meta$condition, label=meta$individual)) + 
  geom_point() + 
  geom_text(check_overlap = TRUE, hjust=0, nudge_x = 0.1) +
  labs(colour='condition') +
  expand_limits(x = c(1, 10))
dev.off()

# K-means clustering on the PCA data *Males more common in study*
pc <- pca$x
km.res <- kmeans(pc[,1:2], centers = 2, nstart = 50)
female <- unique(which.min(table(km.res$cluster)))
sex.list <- ifelse(km.res$cluster == female, 'F', 'M')
# Correcting annotation
sex.list[1] <- 'M'

# Add sex to metadata
pbmc$sex <- sex.list[pbmc$individual]

# Output all cells
saveRDS(pbmc, 'pbmc.RDS')
# Subset for females and output 
pbmc.female <- subset(pbmc, sex == 'F')

pdf('DimPlot.female.pdf')
DimPlot(pbmc.female, label = TRUE, reduction='ref.umap', repel=T) + NoLegend()
dev.off()

DefaultAssay(pbmc) <- 'SCT'
saveRDS(pbmc.female, 'pbmc.female.RDS')