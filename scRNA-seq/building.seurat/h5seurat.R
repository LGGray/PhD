library(Seurat)
library(SeuratDisk)

Convert("local.h5ad", dest="h5seurat", overwrite=T)

pbmc <- LoadH5Seurat("local.h5seurat")

import anndata
import loompy
adata = anndata.read("datasets/lupus.Chun/local.h5ad")
adata.write_loom("datasets/lupus.Chun/local.h5ad.loom", write_obsm_varm=False)

