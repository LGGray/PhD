# load packages
library(Seurat)
library(MAST)
library(scater)
library(data.table)


# read in cell_types.RDS
data <- readRDS("/directflow/SCCGGroupShare/projects/sarba2/data/onek1k/samples/cell_type.RDS")
# read in metadata file
load("/directflow/SCCGGroupShare/projects/sarba2/data/onek1k/samples/metadata.Rdata")
# read in metadata2 file
load("/directflow/SCCGGroupShare/projects/sarba2/data/onek1k/samples/metadata2.Rdata")

data@meta.data = metadata

ids <- c("693_694", "705_706", "733_734", "773_774", "404_405", "866_867",
         "548_549", "590_591", "224_225", "222_223", "692_693",	"687_688",	"708_709",	"983_984",	"984_985",	"749_750",	
         "772_773",	"767_768",	"494_495",	"1066_1067",	"342_343",	"559_560",	
         "393_394",	"602_603", "584_585",	"644_645",	"585_586",	"10_10", "11_11",	
         "9_9",	"220_221",	"221_222",	"215_216")

RA <- c("693_694", "705_706", "733_734", "773_774", "404_405", "866_867",
        "548_549", "590_591", "224_225", "222_223")

control <- c("692_693",	"687_688",	"708_709",	"983_984",	"984_985",	"749_750",	
               "772_773",	"767_768",	"494_495",	"1066_1067",	"342_343",	"559_560",	
               "393_394",	"602_603", "584_585",	"644_645",	"585_586",	"10_10", "11_11",	
               "9_9",	"220_221",	"221_222",	"215_216")

metadata_filtered <- subset(metadata2, individual %in% ids)

pbmc <- subset(data, individual %in% metadata_filtered$individual)

rm(data)

# add disease state to data
for (id in RA){
  for (i in grep(id, pbmc$individual)){
    pbmc@meta.data[i,20] <- "Y"
  }
}

for (id in control){
  for (i in grep(id, pbmc$individual)){
    pbmc@meta.data[i,20] <- "N"
  }
}

colnames(pbmc@meta.data)[20] <- "RA"

print("# of RA")
print(length(grep("Y", unique(paste(pbmc$individual, pbmc$RA, sep="-")))))
print("# of Controls")
print(length(grep("N", unique(paste(pbmc$individual, pbmc$RA, sep="-")))))
print(pbmc)

# Change celltypes to citeseq annotation
Idents(pbmc) <- pbmc$predicted.celltype.l2

# Convert to SingleCellExperiment object
pbmc.sce <- as.SingleCellExperiment(pbmc, assay = "SCT")

# Create a SingleCellAssay object
pbmc.sca = new("SingleCellAssay", pbmc.sce)

# save pbmc file
saveRDS(pbmc, file = "/home/lacgra/datasets/OneK1k/C_vs_RA_DEout/pbmc.RDS")

# Add primerid column to rowData
rowData(pbmc.sca)$primerid <- rownames(rowData(pbmc.sca))

# Differential expression
colData(pbmc.sca)$wellKey <- colnames(pbmc.sca)
cond <- factor(colData(pbmc.sca)$RA)
cond<-relevel(cond, 1)
colData(pbmc.sca)$RA <- cond
zlmCond <- zlm(~RA, pbmc.sca)
print("success")

# Save output
saveRDS(zlmCond, file = "/home/lacgra/datasets/OneK1k/C_vs_RA_DEout/zlmCond.RDS")

# make contrasts
summaryCond_RA <- summary(zlmCond, doLRT='RAY', logFC=TRUE)
saveRDS(summaryCond_RA, file = "/home/lacgra/datasets/OneK1k/C_vs_RA_DEout/summaryCond_RA.RDS")

summaryDt_RA <- summaryCond_RA$datatable
write.table(summaryDt_RA, file = "/home/lacgra/datasets/OneK1k/C_vs_RA_DEout/summaryDt_RA.txt",
            sep="\t", row.names = FALSE)

# Replace column name for convenience
colnames(summaryDt_RA)[4] <- "pvalue"

# Create fcHurdleSig
fcHurdle <- merge(subset(summaryDt_RA, contrast=='RAY' & component=='H')[,c(1,4)], #hurdle P values
                  subset(summaryDt_RA, contrast=='RAY' & component=='logFC')[,c(1, 7, 5, 6)], by='primerid')
fcHurdle$fdr <- p.adjust(fcHurdle$pvalue, 'fdr')
fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(pbmc.sca)), by='primerid')
setorder(fcHurdleSig, pvalue)

write.table(fcHurdleSig, file = paste0("/home/lacgra/datasets/OneK1k/C_vs_RA_DEout/fcHurdleSig.txt"), 
            sep = "\t", row.names = F)






