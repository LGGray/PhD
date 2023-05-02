library(iasva)
library(sva)
library(irlba)
library(Seurat)
library(SummarizedExperiment)

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/pSS_GSE157278')

pbmc <- readRDS('pbmc.female.RDS')

# Calculate geometric library size
geo_lib_size <- colSums(log(pbmc_subset@assays$RNA@data +1))

# Run IA-SVA
set.seed(100)
individual <- pbmc$individual
mod <- model.matrix(~individual + geo_lib_size)
# create a SummarizedExperiment class
sce <- SummarizedExperiment(assay=as.matrix(pbmc@assays$RNA@data))
iasva.res <- iasva(sce, mod[, -1], num.sv = 5)

saveRDS(iasva.res, 'iasva.res.RDS')
