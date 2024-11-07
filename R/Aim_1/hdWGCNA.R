# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# load the Zhou et al snRNA-seq dataset
seurat_obj <- readRDS('Zhou_2020.rds')

# Plot clustering
p <- DimPlot(seurat_obj, group.by='cell_type', label=TRUE) +
   umap_theme() + ggtitle('Zhou et al Control Cortex') + NoLegend()

pdf('Zhou_2020_Clustering.pdf')
print(p)
dev.off()

# Set up the Seurat object for hdWGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)

# Construct metacells - group by cell type and sample
# use k values between 20-75
# Filter out low frequency cell types with min_cells
# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type", "Sample"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cell_type' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

# Set expression matrix for one cell type
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "INH", # the name of the group of interest in the group.by column
  group.by='cell_type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# # Set expression matrix for multiple cell types
# seurat_obj <- SetDatExpr(
#   seurat_obj,
#   group_name = c("INH", "EX"),
#   group.by='cell_type'
# )

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
pdf('Zhou_2020_SoftPowers.pdf')
wrap_plots(plot_list, ncol=2)
dev.off()

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'INH', # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE, # overwrite the TOM if it already exists
)

pdf('Zhou_2020_dendrogram.pdf')
PlotDendrogram(seurat_obj, main='INH hdWGCNA Dendrogram')
dev.off()

# ScaleData before module detection
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
 seurat_obj,
 group.by.vars="Sample"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

save(hMEs, file='Zhou_2020_hMEs.rds')

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type', group_name = 'INH'
)

# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
save(modules, file='Zhou_2020_modules.rds')

# get top N hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
write.csv(hub_df, 'Zhou_2020_hub_genes.csv')

# save the hdWGCNA object
saveRDS(seurat_obj, file='hdWGCNA_object.rds')



