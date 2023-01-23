library(GEOquery)
library(Seurat)

GEO <- "GSE109449"

dir.create(paste0("~/datasets/", GEO))
setwd(paste0("~/datasets/", GEO))

geoFile <- getGEOSuppFiles(GEO, makeDirectory=FALSE)
gzFile <- grep(GEO, basename(rownames(geoFile)), value=TRUE)
txtFile <- gsub(".gz", "", gzFile)
for (file in txtFile){
  gunzip(paste0(file, ".gz"), destname=file, remove=TRUE, overwrite=T)
}

counts <- dir()[grep("gene_counts", dir())]
metadata <- dir()[grep("metadata", dir())]
expr <- read.table(counts)
metadata <- read.table(metadata, header=T)

pbmc <- CreateSeuratObject(expr)
pbmc@meta.data <- metadata

saveRDS(pbmc, "pbmc.RDS")

# counts <- read.delim("datasets/SDY998/celseq_matrix_ru1_reads.tsv.725586", row.names=1, header=T)
# counts[is.na(counts)] <- 0
# meta <- read.delim("datasets/SDY998/celseq_meta.tsv.725591", header=T)
# meta <- subset(meta, type != "Empty")
# counts[1:3,1:3]
# 
# pbmc <- CreateSeuratObject(counts)
# pbmc@meta.data <- meta


# pbmc <- NormalizeData(object = pbmc)
# pbmc <- FindVariableFeatures(object = pbmc)
# pbmc <- ScaleData(object = pbmc)
# pbmc <- RunPCA(object = pbmc)
# pbmc <- FindNeighbors(object = pbmc)
# pbmc <- FindClusters(object = pbmc)
# pbmc <- RunTSNE(object = pbmc)
# pdf("plots/SDY998/dimplot.pdf")
# DimPlot(object = pbmc, reduction = "tsne")
# dev.off()


