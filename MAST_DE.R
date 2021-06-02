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

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

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

# Filter out NA values in metadata
tokeep <- rownames(metadata[complete.cases(metadata), ])
data_subset <- data_subset[,colnames(data_subset) %in% tokeep]

# Change pool_X to numeric
data_subset$pool <- as.numeric(sub("pool_", "", data_subset$pool))

# Change celltypes to citeseq annotation
Idents(data_subset) <- data_subset$predicted.celltype.l2

# Select a number of pools to analyze
pbmc <- subset(data_subset, pool %in% args[1]:args[1]+10)
print("# of Males")
print(length(grep("-1", unique(paste(pbmc$individual, pbmc$sex, sep="-")))))
print("# of Females")
print(length(grep("-2", unique(paste(pbmc$individual, pbmc$sex, sep="-")))))
print(pbmc)

# Convert to SingleCellExperiment object
pbmc.sce <- as.SingleCellExperiment(pbmc, assay = "SCT")

# Create a SingleCellAssay object
pbmc.sca = new("SingleCellAssay", pbmc.sce)

# save file
saveRDS(zlmCond, file = paste0("/home/lacgra/datasets/OneK1k/M_vs_F_DEout/zlmCond_RDS/zlmCond", args[1], "-", args[1]+10, ".RDS")
saveRDS(pbmc, file = paste0("/home/lacgra/datasets/OneK1k/M_vs_F_DEout/pbmc_RDS/pbmc", args[1], "-", args[1]+10, ".RDS")

# Differential expression 
options(mc.cores = 4)
zlmCond <- zlm(~sex + predicted.celltype.l2 + (1 | pool) + (1 | age), pbmc.sca, method= 'bayesglm', parallel = TRUE)
summaryCond_sex <- summary(zlmCond, doLRT='sex')
summaryCond_celltype <- summary(zlmCond, doLRT='predicted.celltype.l2')
summaryCond_all <- summary(zlmCond, doLRT=TRUE)

summaryDt_sex <- summaryCond_sex$datatable
summaryDt_celltype <- summaryCond_celltype$datatable
summaryDt_all <- summaryCond_all$datatable

write.table(summaryDt_sex, file=paste0("/home/lacgra/datasets/OneK1k/MAST_out/summaryDt_sex", args[1], "-", args[1]+10, ".txt"), sep="\t", row.names = FALSE)
write.table(summaryDt_celltype, file=paste0("/home/lacgra/datasets/OneK1k/MAST_out/summaryDt_celltype", args[1], "-", args[1]+10, ".txt"), sep="\t", row.names = FALSE)
write.table(summaryDT_all, file=paste0("/home/lacgra/datasets/OneK1k/MAST_out/summaryDt_all", args[1], "-", args[1]+10, ".txt"), sep="\t", row.names = FALSE

