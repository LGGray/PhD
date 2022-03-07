suppressPackageStartupMessages({
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(KernSmooth)
  library(BiocParallel)
  library(ggplot2)
  library(data.table)
  library(grid)
  library(ComplexHeatmap)
  library(Seurat)
  library(doParallel)
  library(doRNG)
  library(viridis)
  library(scFunctions)
  library(reshape2)
  library(ggridges)
  library(tidyverse)
})

options(width=200)

# setwd("~/datasets/OneK1k/co-express/SCENIC/all/")

# # # Read in data and extract expression matrix and cluster information
# pbmc <- readRDS("/directflow/SCCGGroupShare/projects/lacgra/seurat.object/pbmc.all.RDS")
# pbmc <- subset(pbmc, features = VariableFeatures(pbmc))
# exprMat <- GetAssayData(pbmc, slot="counts")
# exprMat <- as.matrix(exprMat)
# saveRDS(exprMat, "~/datasets/OneK1k/co-express/SCENIC/all/int/exprMat.Rds")
# cellInfo <- data.frame(seuratCluster=Idents(pbmc))
# saveRDS(cellInfo, "~/datasets/OneK1k/co-express/SCENIC/all/int/cellInfo.Rds")
# #
# # ### Initialize settings
# library(SCENIC)
# scenicOptions <- initializeScenic(org="hgnc", dbDir="~/genome.files/cisTarget_databases", nCores=8)
# saveRDS(scenicOptions, file="int/scenicOptions.Rds")
# #
# # ### Co-expression network
# genesKept <- geneFiltering(exprMat, scenicOptions)
# exprMat_filtered <- exprMat[genesKept, ]
# runCorrelation(exprMat_filtered, scenicOptions)
# exprMat_filtered_log <- log2(exprMat_filtered+1)
# #runGenie3(exprMat_filtered_log, scenicOptions)
# exportsForArboreto(exprMat_filtered_log, scenicOptions, dir = "int")

### load the arboreto output and save as GENIE3 output
setwd("~/datasets/OneK1k/co-express/SCENIC/all/")
linkList <- read.delim("int/arboreto.tsv", header=F)
colnames(linkList) <- c("TF", "Target", "weight")
saveRDS(linkList, "int/1.4_GENIE3_linkList.Rds")

#--------------------------------------------------------------------------

# If starting a new session following arboreto
setwd("~/datasets/OneK1k/co-express/SCENIC/all/")
exprMat <- readRDS("int/exprMat.Rds")
cellInfo <- readRDS("int/cellInfo.Rds")
genesKept <- readRDS("int/1.1_genesKept.Rds")
exprMat_filtered <- read.delim("int/1.1_exprMatrix_filtered_t.txt")
exprMat_filtered <- log2(exprMat_filtered+1)

scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
#scenicOptions@settings$dbs <- scenicOptions@settings$dbs[("10kb")]
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# Optional: Binarize activity
# scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
# tsneAUC(scenicOptions, aucType="AUC") # choose settings

# Export:
saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="XBP1"]
viewMotifs(tableSubset)

# output/Step2_regulonTargetsInfo.tsv in detail:
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="STAT6" & highConfAnnot==TRUE]
viewMotifs(tableSubset)

# Cell-type specific regulators (RSS):
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "seuratCluster"])
rssPlot <- plotRSS(rss)
pdf("output/RSS.plot.pdf")
rssPlot$plot
dev.off()

# #-------------------------------------------------------------------------
# X-linked and escape Regulon targets
load("~/datasets/OneK1k/X_escape/chrX.Rdata")
load("~/datasets/OneK1k/X_escape/escapees.Rdata")
regulonTargetsInfo <- subset(regulonTargetsInfo, highConfAnnot == T)
X.targets <- subset(regulonTargetsInfo, gene %in% rownames(chrX))
escape.targets <- subset(regulonTargetsInfo, gene %in% rownames(escape))
# Select TFs which target XCI escapees
features <- unique(escape.targets$TF)
# 
# #-------------------------------------------------
# Load in seurat object
pbmc <- readRDS("/directflow/SCCGGroupShare/projects/lacgra/seurat.object/pbmc.SCENIC.RDS")
# Add data to Seurat object
all.equal(colnames(pbmc), colnames(regulonAUC))
AUCmat <- AUCell::getAUC(regulonAUC)
# Remove low confidence TFs from AUC, AUCbinary and rss
AUCmat <- AUCmat[!grepl("_extended", rownames(AUCmat)),]
rownames(AUCmat) <- gsub('[()]|\\d+g| ', '',rownames(AUCmat))
pbmc[['AUC']] <- CreateAssayObject(data = AUCmat)
# #----
# binaryMat <- as.matrix(binaryMat)
# binaryMat <- binaryMat[!grepl("_extended", rownames(binaryMat)),]
# rownames(binaryMat) <- gsub('[()]|\\d+g| ', '',rownames(binaryMat))
# pbmc[['AUCBinary']] <- CreateAssayObject(data = binaryMat)
# #----
rss <- rss[!grepl("_extended", rownames(rss)),]
rownames(rss) <- gsub('[()]|\\d+g| ', '',rownames(rss))
pbmc[['rss']] <- CreateAssayObject(data = rss)

saveRDS(pbmc, "/directflow/SCCGGroupShare/projects/lacgra/seurat.object/pbmc.SCENIC.RDS")
# 
# Visualizing data
setwd("~/plots/OneK1K/co-express/SCENIC/all")
#---------------------------------
DefaultAssay(pbmc) <- 'AUC'
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(pbmc)
pdf("AUC.heatmap.celltype.pdf", width = 20, height = 14)
DoHeatmap(pbmc, features=features, slot = 'scale.data', size=3) +
  scale_fill_viridis()
dev.off()
pdf("NFKB1.ridgeplot.celltype.AUC.pdf")
RidgePlot(pbmc, feature = "NFKB1", split.by="RA", same.y.lims=T)
dev.off()
pdf("DEG.dotplot.celltype.AUC.pdf")
DotPlot(pbmc, features= c("JUN", "STAT1", "FOS", "GTF2F1", "RUNX3", "NR3C1"), split.by="RA") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Cell type") +
  xlab("DE Regulon") +
  ggtitle("Differentially Expressed regulon AUC")
dev.off()
# 
# 
# #-----------------------
# DefaultAssay(pbmc) <- 'AUCBinary'
# pdf("AUCBinary.heatmap.condition.pdf", width = 20, height = 14)
# DoHeatmap(pbmc, features=features, group.by = 'RA', slot = 'data', size=3) +
#   scale_fill_viridis()
# dev.off()
# pdf("NFKB1.ridgeplot.AUCBinary.pdf")
# RidgePlot(pbmc, feature = "NFKB1")
# dev.off()
# 
# #------------------------------
# # Using scFunctions
# # Regulon Specificity Score
# setwd("~/plots/OneK1K/co-express/SCENIC/all/")
# regulonAUC <- readRDS("~/datasets/OneK1k/co-express/SCENIC/all/int/3.4_regulonAUC.Rds")
# kmeans_thresholds <- auc_thresh_kmeans(regulonAUC)
# binary_regulons <- binarize_regulons(regulonAUC,kmeans_thresholds)
# joined_bin_reg <- binary_regulons %>%
#   purrr::reduce(left_join,by="cells") %>%
#   tibble::column_to_rownames("cells")
# 
# binary_regulons_trans <- as.matrix(t(joined_bin_reg))
# rrs_df <- melt(readRDS("~/datasets/OneK1k/co-express/SCENIC/all/output/rss.RDS"))
# colnames(rrs_df) <- c("regulon", "cell_type", "RSS")
# pdf("RRS_rankings.all.pdf")
# plot_rrs_ranking(rrs_df,
#                  "all",
#                  ggrepel_force = 1,
#                  ggrepel_point_padding = 0.2,
#                  top_genes = 4,
#                  plot_extended = FALSE)
# dev.off()
# rrs_df_deg <- subset(rrs_df, regulon %in% features)
# pdf("RRS_rankings.DEG.pdf")
# plot_rrs_ranking(rrs_df_deg,
#                  "all",
#                  ggrepel_force = 1,
#                  ggrepel_point_padding = 0.2,
#                  top_genes = 6,
#                  plot_extended = FALSE)
# dev.off()
# 
# pdf("RSS.ridgeplot.pdf")
# rrs_df_nona <- subset(rrs_df,RSS > 0)
# ggplot(rrs_df_nona,aes(RSS,cell_type, fill = cell_type)) +
#   geom_density_ridges(scale = 5, alpha = 0.75) +
#   geom_vline(xintercept = 0.1) +
#   theme(legend.position = "none")
# dev.off()
# pdf("RSS.ridgeplot.DEG.pdf")
# rrs_df_nona <- subset(rrs_df_deg,RSS > 0)
# ggplot(rrs_df_nona,aes(RSS,cell_type, fill = cell_type)) +
#   geom_density_ridges(scale = 5, alpha = 0.75) +
#   geom_vline(xintercept = 0.1) +
#   theme(legend.position = "none")
# dev.off()
# 
# 
# rrs_df_wide <- rrs_df %>%
#   spread(cell_type,RSS)
# 
# rownames(rrs_df_wide) <- rrs_df_wide$regulon
# rrs_df_wide <- rrs_df_wide[,2:ncol(rrs_df_wide)]
# 
# ## Subset all regulons that don't have at least an RSS of 0.7 for one cell type
# rrs_df_specific <- subset(rrs_df, RSS < 0.3)
# rrs_df_specific
# 
# pdf("RSS.heatmap.pdf", height=15)
# ggplot(data=rrs_df_specific, aes(x = cell_type, y=regulon, fill = RSS)) + geom_tile() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
# dev.off()
# 
# # ----------------------------
# # Connection Specificity Index
# regulons_csi <- calculate_csi(regulonAUC,
#                               calc_extended = FALSE)
# pdf("CSI.modules.pdf")
# plot_csi_modules(regulons_csi,
#                  nclust = 10,
#                  font_size_regulons = 8)
# dev.off()
# 
# regulons_csi$regulon_1 <- gsub('[()]|\\d+g| ', '', regulons_csi$regulon_1)
# regulons_csi$regulon_2 <- gsub('[()]|\\d+g| ', '', regulons_csi$regulon_2)
# regulons_csi_deg <- subset(regulons_csi, regulon_1 %in% features &
#                              regulon_2 %in% features)
# pdf("CSI.modules.DEG.pdf")
# ggplot(regulons_csi_deg, aes(x=regulon_1, y=regulon_2, fill = CSI)) +
#   geom_tile() +
#   scale_fill_viridis() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=6),
#         axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=6)) +
#   ggtitle("CSI Differentially Expressed regulons") +
#   labs(fill="CSI")
# 
# 
# # Export clusters
# csi_csi_wide <- regulons_csi %>%
#   spread(regulon_2,CSI)
# 
# future_rownames <- csi_csi_wide$regulon_1
# csi_csi_wide <- as.matrix(csi_csi_wide[,2:ncol(csi_csi_wide)])
# rownames(csi_csi_wide) <- future_rownames
# 
# regulons_hclust <- hclust(dist(csi_csi_wide,method = "euclidean"))
# 
# clusters <- cutree(regulons_hclust,k= 10)
# clusters_df <- data.frame("regulon" = names(clusters),
#                           "csi_cluster" = clusters)
# # Check how many regulons are in each cluster
# clusters_df_stats <- clusters_df %>%
#   group_by(csi_cluster) %>%
#   mutate("regulon" = as.character(regulon)) %>%
#   tally()
# pdf("CSI.clusters.pdf")
# ggplot(clusters_df_stats,aes(as.factor(csi_cluster),n,fill=as.factor(csi_cluster))) +
#   geom_bar(color= "black",stat="identity") +
#   theme(legend.position="none") +
#   scale_fill_brewer(palette = "Set3") +
#   labs(x = "HC clusters",
#        y = "# Regulons")
# dev.off()
# ## Check average regulon size per cluster
# clusters_df_regsizes <- clusters_df %>%
#   separate(regulon, into = c("regulon_name","regulon_size"), sep=" ") %>%
#   mutate("regulon_size" = gsub("\\(","",regulon_size)) %>%
#   mutate("regulon_size" = gsub("\\g)","",regulon_size)) %>%
#   mutate("regulon_size" = as.numeric(regulon_size))
# pdf("CSI.regulon.cluster.pdf")
# ggplot(clusters_df_regsizes,aes(log10(regulon_size),as.factor(csi_cluster),fill=as.factor(csi_cluster))) +
#   geom_density_ridges() +
#   scale_fill_brewer(palette = "Set3") +
#   theme(legend.position = "none") +
#   labs(x = "log10 regulon size", y = "csi cluster")
# dev.off()
