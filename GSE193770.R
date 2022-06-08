# Script for building Seurat object from GSE193770 data

library(Seurat)
library(ddqcR)
library(SeuratDisk)

setwd('~/datasets/GSE193770')

files <- list.files(patter='counts.txt', recursive = T, full.names = T)
data <- lapply(files, function(x) read.delim(x))
names(data) <- lapply(files, function(x){
  unlist(strsplit(basename(x), '_'))[2]
})

# Simple test for sex
lapply(data, function(x){
  rowSums(subset(x, rownames(x) %in% c('XIST', 'RPS4Y', 'RPS4X', 'DDX3Y')))
})

sex <- c('F', 'F', 'M', 'F', 'F', 'M', 'M', 'F', 'F', 'M')

# Create metadata
metadata <- lapply(seq_along(data), function(x){
  data.frame(cell_id=colnames(data[[x]]), 
                     sample=rep(names(data)[[x]], ncol(data[[x]])),
                     sex=rep(sex[x], ncol(data[[x]])))
})

metadata <- do.call("rbind", metadata)
rownames(metadata) <- metadata$cell_id

# Bind columns across samples
exprMat <- do.call("cbind", data)
# Edit cell names
colnames(exprMat) <- gsub('MS[0-9]+\\.|HC[0-9]+\\.', '', colnames(exprMat))

# Create Seurat object
pbmc <- CreateSeuratObject(counts=exprMat)
meta <- cbind(pbmc@meta.data, metadata)

# Check that cell names are equal
all.equal(gsub('\\.[0-9]', '', colnames(pbmc)), gsub('\\.[0-9]', '', pbmc@meta.data$cell_id))

# Remove obvious bad quality cells
pbmc <- initialQC(pbmc)

# Return dataframe of filtering statistics
pdf('ddqc.plot.pdf')
df.qc <- ddqc.metrics(pbmc)
dev.off()

# Filter out the cells
pbmc <- filterData(pbmc, df.qc)

# SCTransform
pbmc <- SCTransform(pbmc, verbose = FALSE)

# Read in PBMC reference dataset
reference <- LoadH5Seurat("~/azimuth.reference/pbmc_multimodal.h5seurat")
# Find anchors between reference and query
anchors <- FindTransferAnchors(reference = reference,
                               query = pbmc,
                               normalization.method = "SCT",
                               reference.reduction = "spca",
                               dims = 1:50)
# Transfer cell type labels and protein data from reference to query
pbmc <- MapQuery(anchorset = anchors,
                 query = pbmc,
                 reference = reference,
                 refdata = list(celltype.l1 = "celltype.l1",
                                celltype.l2 = "celltype.l2",
                                predicted_ADT = "ADT"),
                 reference.reduction = "spca",
                 reduction.model = "wnn.umap")

Idents(pbmc) <- 'predicted.celltype.l2'

pbmc@meta.data$condition <- gsub('[0-9]+', '', pbmc@meta.data$sample)

# Checking cell counts per condition
MS <- subset(pbmc, condition == 'MS')
HC <- subset(pbmc, condition == 'HC')

ms <- data.frame(MStable(MS$predicted.celltype.l2))
colnames(ms) <- c('celltype', 'freq')
hc <- data.frame(HC=table(HC$predicted.celltype.l2))
colnames(hc) <- c('celltype', 'freq')
plot.data <- merge(ms, hc, by='celltype', all=T)
colnames(plot.data) <- c('celltype', 'MS','HC')
plot.data[is.na(plot.data)] <- 0
plot.data <- reshape2::melt(plot.data)
ggplot(plot.data, aes(x=celltype, y=value, colour=variable)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle=45, hjust=1))

cells <- levels(pbmc)
keep <- c('CD4 CTL', 'CD4 Naive', 'CD4 TCM', 'CD4 TEM', 
          'CD8 Naive', 'CD8 TCM', 'CD8 TEM','dnT', 'gdT', 
          'MAIT', 'NK', 'Treg')
pbmc <- subset(pbmc, predicted.celltype.l2 %in% keep)

pdf('DimPlot.pdf')
DimPlot(pbmc, label = TRUE, reduction='ref.umap', repel=T) + NoLegend()
dev.off()
saveRDS(pbmc, 'pbmc.female.RDS')
