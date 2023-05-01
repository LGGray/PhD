library(iasva)
library(Seurat)
library(SummarizedExperiment)

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/pSS_GSE157278')

pbmc <- readRDS('pbmc.female.RDS')

# Calculate geometric library size
geo_lib_size <- colSums(log(pbmc@assays$RNA@data +1))

# Run IA-SVA
set.seed(100)
individual <- pbmc$individual
mod <- model.matrix(~individual + geo_lib_size)
# create a SummarizedExperiment class
sce <- SummarizedExperiment(assay=as.matrix(pbmc@assays$RNA@data))
iasva.res <- fast_iasva(sce, mod[, -1], num.sv = 5)

saveRDS(iasva.res, 'iasva.res.RDS')
