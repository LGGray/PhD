# Script information ----------------------------------------------------------

# title: Subset the OneK1k data to control individuals
# Author: Lachlan G Gray
# Date: 17/05/2021
# Description: Convert the OneK1k Seurat object to a SingleCellExperiment object. Then perform a differential expression analysis with MAST.
# Filter for 11 pools per job
# Find genes differentally expressed between males and females for each celltype, factoring in age and batch effects.
# Run as an array job

# load packages
library(Seurat)
library(MAST)
library(scater)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

# read in cell_types.RDS
data <- readRDS("/directflow/SCCGGroupShare/projects/sarba2/data/onek1k/samples/cell_type.RDS")
# read in metadata file
load("/directflow/SCCGGroupShare/projects/sarba2/data/onek1k/samples/metadata.Rdata")
# read in metadata2 file
load("/directflow/SCCGGroupShare/projects/sarba2/data/onek1k/samples/metadata2.Rdata")

data@meta.data = metadata

# print object
print(data)

# filter metadata file
metadata_filtered <- metadata2[metadata2$X08_Autoimmune_Disease == "N" & 
                                 metadata2$X11_Rheumatoid_arthritis == "N" & 
                                 metadata2$X17_Cancer == "N",]
metadata_filtered <- metadata_filtered[ !grepl("NA", rownames(metadata_filtered)), ]

# subset data based on filtered metadata
data_subset <- subset(data, individual %in% metadata_filtered$individual)

rm(data)

# Filter out NA values in metadata
tokeep <- rownames(metadata[complete.cases(metadata), ])
data_subset <- data_subset[,colnames(data_subset) %in% tokeep]

# Change pool_X to numeric
data_subset$pool <- as.numeric(sub("pool_", "", data_subset$pool))

# Change celltypes to citeseq annotation
Idents(data_subset) <- data_subset$predicted.celltype.l2

pbmc <- subset(data_subset, pool %in% args[1])
print("# of Males")
print(length(grep("-1", unique(paste(pbmc$individual, pbmc$sex, sep="-")))))
print("# of Females")
print(length(grep("-2", unique(paste(pbmc$individual, pbmc$sex, sep="-")))))
print(pbmc)

rm(data_subset)

# Convert to SingleCellExperiment object
pbmc.sce <- as.SingleCellExperiment(pbmc, assay = "SCT")

# Create a SingleCellAssay object
pbmc.sca = new("SingleCellAssay", pbmc.sce)

# save pbmc file
saveRDS(pbmc, file = paste0("/home/lacgra/datasets/OneK1k/M_vs_F_DEout/pbmc_RDS/pbmc", args[1], ".RDS"))

# Add primerid column to rowData
rowData(pbmc.sca)$primerid <- rownames(rowData(pbmc.sca))

# Differential expression 
colData(pbmc.sca)$wellKey <- colnames(pbmc.sca)
cond <- factor(colData(pbmc.sca)$sex)
cond<-relevel(cond, 1)
colData(pbmc.sca)$sex <- cond
zlmCond <- zlm(~sex + predicted.celltype.l2, pbmc.sca)
print("success")

# save file
saveRDS(zlmCond, file = paste0("/home/lacgra/datasets/OneK1k/M_vs_F_DEout/zlmCond_RDS/zlmCond", args[1], ".RDS"))

# make contrasts
summaryCond_sex <- summary(zlmCond, doLRT='sex2', logFC=TRUE)
#summaryCond_all <- summary(zlmCond, doLRT=TRUE, logFC=TRUE)

saveRDS(summaryCond_sex, file = paste0("/home/lacgra/datasets/OneK1k/M_vs_F_DEout/summaryCond_RDS/summaryCond_sex", args[1], ".RDS"))
#saveRDS(summaryCond_all, file = paste0("/home/lacgra/datasets/OneK1k/M_vs_F_DEout/summaryCond_RDS/summaryCond_all", args[1], ".RDS"))

summaryDt_sex <- summaryCond_sex$datatable
#summaryDt_all <- summaryCond_all$datatable
#colnames(summaryDt_all)[4] <- "logFC"

write.table(summaryDt_sex, file=paste0("/home/lacgra/datasets/OneK1k/M_vs_F_DEout/MAST_out/summaryDt_sex", args[1], ".txt"), sep="\t", row.names = FALSE)
#write.table(summaryDt_all, file=paste0("/home/lacgra/datasets/OneK1k/M_vs_F_DEout/MAST_out/summaryDt_all", args[1], ".txt"), sep="\t", row.names = FALSE)

# Replace column name for convenience
colnames(summaryDt_sex)[4] <- "pvalue"

# Create fcHurdleSig
fcHurdle <- merge(subset(summaryDt_sex, contrast=='sex2' & component=='H')[,c(1,4)], #hurdle P values
                  subset(summaryDt_sex, contrast=='sex2' & component=='logFC')[,c(1, 7, 5, 6)], by='primerid')
fcHurdle$fdr <- p.adjust(fcHurdle$pvalue, 'fdr')
fcHurdleSig <- merge(subset(fcHurdle, fdr < 0.05), as.data.table(mcols(pbmc.sca)), by='primerid')
setorder(fcHurdleSig, fdr)

write.table(fcHurdleSig, file = paste0("fcHurdleSig/pool_", args[1], ".txt"), 
            sep = "\t", row.names = F)

# Visualize top 50 DEG
filename <- paste0('pool_', args[1], ".pdf")
pdf(file.path("~/plots/DE_genes", filename), height = 14, width = 14)
entrez_to_plot <- fcHurdleSig[1:50,1]
flat_dat <- as(pbmc.sca[entrez_to_plot,], 'data.table')
ggbase <- ggplot(flat_dat, aes(x=sex, y=logcounts,color=sex)) + geom_jitter()+facet_wrap(~primerid, scale='free_y')+ggtitle("DE Genes")
ggbase+geom_violin()
dev.off()

