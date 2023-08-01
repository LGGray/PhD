library(Seurat)
library(SeuratDisk)
library(loomR)

# if pbmc.loom exists then end the script
if (file.exists("pbmc.loom")) {
  print("pbmc.loom already exists")
  q()
} else {
  pbmc <- readRDS('pbmc.female.RDS')
  pbmc <- subset(pbmc, features = VariableFeatures(pbmc))

  expr <-GetAssayData(pbmc, slot = "counts")

  # Create a new loom file
  lfile <- create(filename = "pbmc.loom", data = expr)

  # Add the row and column attributes
  lfile$addAttr("CellID", colnames(expr), axis = 1)
  lfile$addAttr("Gene", rownames(expr), axis = 0)

  # Close the connection to the loom file
  lfile$disconnect()
}

TF <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/SCENIC/allTFs_hg38.txt', header=F)
print('# of matches '+ sum(rownames(pbmc) %in% TF$V1))

