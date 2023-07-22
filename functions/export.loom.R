library(Seurat)
library(SeuratDisk)
library(loomR)

# if pbmc.loom exists then end the script
if (file.exists("pbmc.loom")) {
  print("pbmc.loom already exists")
  q()
} else {
  pbmc <- readRDS('pbmc.female.RDS')

  expr <-GetAssayData(pbmc, slot = "counts")

  # Create a new loom file
  lfile <- create(filename = "pbmc.loom", data = expr)
  # Close the connection to the loom file
  lfile$disconnect()
}



