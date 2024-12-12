# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# gene enrichment packages
library(enrichR)
library(GeneOverlap)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# Create dir
if (!dir.exists('hdWGCNA')) {
  dir.create('hdWGCNA')
}

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# load the Zhou et al snRNA-seq dataset
seurat_obj <- readRDS('pbmc.female.RDS')

# Set up the Seurat object for hdWGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "wgcna" # the name of the hdWGCNA experiment
)

# Construct metacells - group by cell type and sample
# use k values between 20-75
# Filter out low frequency cell types with min_cells
# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cellTypist", "individual"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 30, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cellTypist' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

metacell_obj <- GetMetacellObject(seurat_obj)
celltypes <- levels(metacell_obj)

# Set expression matrix for one cell type
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "Tcm/Naive helper T cells", # the name of the group of interest in the group.by column
  group.by='cellTypist', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'decontXcounts', # using RNA assay
  slot = 'data' # using normalized data
)

# # Set expression matrix for multiple cell types
# seurat_obj <- SetDatExpr(
#   seurat_obj,
#   group_name = celltypes,
#   group.by='cellTypist'
# )

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
pdf('hdWGCNA/SoftPowers.pdf')
wrap_plots(plot_list, ncol=2)
dev.off()

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'wgcna', # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE, # overwrite the TOM if it already exists
)

pdf('hdWGCNA/dendrogram.pdf')
PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
dev.off()

TOM <- GetTOM(seurat_obj)
save(TOM, file='hdWGCNA/TOM.rds')

# ScaleData before module detection
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
 seurat_obj,
 group.by.vars="individual"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)
save(hMEs, file='hdWGCNA/hMEs.rds')

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cellTypist', group_name = 'Tcm/Naive helper T cells'
)

# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
save(modules, file='hdWGCNA/modules.rds')

# get top N hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
write.csv(hub_df, 'hdWGCNA/hub_genes.csv')

# save the hdWGCNA object
saveRDS(seurat_obj, file='hdWGCNA_object.rds')


### Differential module eigengene (DME) analysis ###

# seurat_obj <- readRDS('hdWGCNA_object.rds')

group1 <- seurat_obj@meta.data %>% subset(cellTypist == 'Tcm/Naive helper T cells' & condition == 'control') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(cellTypist == 'Tcm/Naive helper T cells' & condition == 'disease') %>% rownames

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group2,
  barcodes2 = group1,
  test.use='wilcox',
  wgcna_name='wgcna'
)

pdf('hdWGCNA/DME_volcano.pdf')
PlotDMEsVolcano(
  seurat_obj,
  DMEs,
  wgcna_name = 'wgcna'
)
dev.off()

dbs <- 'GO_Biological_Process_2021'

seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

enrich_df <- GetEnrichrTable(seurat_obj)

subset(enrich_df, module == 'turquoise' & Adjusted.P.value < 0.05 & Combined.Score > 200)$Term

test <- subset(enrich_df, module == 'brown' & Adjusted.P.value < 0.05)$Term

# Select Term if contains chrX genes
test$nchrX <- lapply(test$Genes, function(x) sum(strsplit(x, ';')[[1]] %in% chrX))
subset(test, nchrX > 0)

# top 10 enriched terms ranked by Combined.Score
test %>% arrange(desc(Combined.Score)) %>% head(10)

