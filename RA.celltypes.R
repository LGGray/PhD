library(Seurat)
library(MAST)
library(data.table)

# Read in object
pbmc <- readRDS("/home/lacgra/datasets/OneK1k/C_vs_RA_DEout/pbmc.RDS")

# Start loop
for (cell in levels(pbmc)){
  
  print(cell)
  pbmc.cell <- subset(pbmc, predicted.celltype.l2 %in% cell)
  
  # Convert to SingleCellExperiment object
  pbmc.sce <- as.SingleCellExperiment(pbmc.cell, assay = "SCT")
  
  # Create a SingleCellAssay object
  pbmc.sca = new("SingleCellAssay", pbmc.sce)
  
  # Add primerid column to rowData
  rowData(pbmc.sca)$primerid <- rownames(rowData(pbmc.sca))
  
  # Differential expression
  colData(pbmc.sca)$wellKey <- colnames(pbmc.sca)
  cond <- factor(colData(pbmc.sca)$RA)
  cond<-relevel(cond, 1)
  colData(pbmc.sca)$RA <- cond
  zlmCond <- zlm(~RA, pbmc.sca)
  print("success")
  
  # make contrasts
  summaryCond_RA <- summary(zlmCond, doLRT='RAY', logFC=TRUE)
  print("success")
  
  summaryDt_RA <- summaryCond_RA$datatable
  
  # Replace column name for convenience
  colnames(summaryDt_RA)[4] <- "pvalue"
  
  # Create fcHurdleSig
  fcHurdle <- merge(subset(summaryDt_RA, contrast=='RAY' & component=='H')[,c(1,4)], #hurdle P values
                    subset(summaryDt_RA, contrast=='RAY' & component=='logFC')[,c(1, 7, 5, 6)], by='primerid')
  fcHurdle$fdr <- p.adjust(fcHurdle$pvalue, 'fdr')
  fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(pbmc.sca)), by='primerid')
  setorder(fcHurdleSig, pvalue)
  
  # save output
  write.table(fcHurdleSig, file = paste0("/home/lacgra/datasets/OneK1k/C_vs_RA_DEout/celltype/", cell, ".txt"), 
              sep = "\t", row.names = F)
}
  
