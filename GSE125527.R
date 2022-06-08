library(Seurat)
library(ddqcR)
library(SeuratDisk)

setwd('datasets/GSE125527')

counts <- read.csv("GSE125527_UMI_cell_table_sparse.csv", header=T)
features <- read.csv("GSE125527_gene_id_rownames.csv", header=F)
cell.names <- read.csv("GSE125527_cell_id_colnames.csv", header=F)
metadata <- read.csv("GSE125527_cell_metadata.csv", row.names = 1, header=T)

for (i in 1:nrow(counts)){
  counts.matrix[counts[i,1], counts[i,2]] <- counts[i,3]
}

rownames(counts.matrix) <- features$V1
colnames(counts.matrix) <- cell.names$V1

pbmc <- CreateSeuratObject(counts.matrix)
pbmc@meta.data <- cbind(pbmc@meta.data, metadata)

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

Idents(pbmc) <- 'predicted.celltype.l2'

AverageExpression(pbmc, features = c("XIST", "DDX3Y", "RPS4Y"), group.by = 'patient_assignment')

mydict <- c("C1"="M", "C2"="M", "C3"="F", "C4"="M", "C5"="M", "C6"="M", "C7"="M", "C8"="F",
            "C9"="F", "U1"="F", "U2"="F", "U3"="M", "U4"="M", "U5"="F", "U6"="F", "U7"="M")
pbmc$sex <- mydict[pbmc$patient_assignment]

pbmc.female <- subset(pbmc, sex == "F")


keep <- names(which(table(pbmc.female$predicted.celltype.l2) > 100))
pbmc.female <- subset(pbmc.female, predicted.celltype.l2 %in% keep)

pdf('DimPlot.pdf')
DimPlot(pbmc.female, label = TRUE, reduction='ref.umap', repel=T) + NoLegend()
dev.off()

saveRDS(pbmc.female, 'pbmc.female.RDS')

keep <- names(which(table(pbmc$predicted.celltype.l2) > 100))
pbmc <- subset(pbmc, predicted.celltype.l2 %in% keep)
saveRDS(pbmc, 'pbmc.RDS')


