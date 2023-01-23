library(Seurat)
library(nebula)
library(ggplot2)
library(ggrepel)

# read in file
pbmc <- readRDS("~/datasets/OneK1k/C_vs_RA_DEout/pbmc.RDS")

# Iterate over each cell type
for (cell in levels(pbmc)){
  print(cell)
  # subset object by cell type
  pbmc.cell <- subset(pbmc, predicted.celltype.l2 %in% cell)
  
  # Extracting data from object
  count <- pbmc.cell@assays$SCT@counts
  sid <- pbmc.cell@meta.data$individual
  cc <- pbmc.cell@meta.data$RA
  
  pred <- data.frame(cc)
  df = model.matrix(~cc, data = pred)
  # group data for nebula
  data_g = group_cell(count=count, id=sid,pred=df)
  # Run nebula analysis
  re = nebula(data_g$count, data_g$id, pred=data_g$pred)
  
  saveRDS(re, paste0("~/datasets/OneK1k/C_vs_RA_DEout/nebula/", cell, ".RDS"))
}

setwd("~/datasets/OneK1k/C_vs_RA_DEout/nebula/")
# read in XCI escape data
load("../../X_escape/escapees.Rdata")
xcape_genes <- rownames(escape)

# Reading in nebula output files
files <- list.files(pattern = "\\.RDS$")

# Set ggrepel 
options(ggrepel.max.overlaps = Inf)
# Create empty plot list
plot_list = list()
# Iterate over nebula files and create volcano plot
for (i in 1:length(files)){
  cell <- gsub(".RDS", "", files[i])
  print(cell)
  # Read in file
  re <- readRDS(files[i])
  # Indicate significance
  threshold <- re$summary$p_ccY <= 0.05 
  re$summary$threshold <- threshold
  # Order on p value
  re$summary_ordered <- re$summary[order(re$summary$p_ccY),]
  # Indicate if gene is XCI escaping and DE
  re$summary_ordered$xcape <- ifelse(re$summary_ordered$gene %in% xcape_genes & 
                                       re$summary_ordered$threshold == T, T, F)
  
  options(ggrepel.max.overlaps = Inf)
  p = ggplot(re$summary_ordered) +
    geom_point(aes(x=logFC_ccY, y=-log10(p_ccY), colour=threshold)) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    geom_text_repel(box.padding = 0.75, aes(x=logFC_ccY, y=-log10(p_ccY), label = ifelse(xcape == T , gene, ""))) +
    ggtitle(paste(cell, "Volcano Plot")) +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
  
  plot_list[[i]] = p

}

pdf("~/plots/RA_genes/celltype.VP.pdf")
for (i in 1:length(files)){
  print(plot_list[[i]])
}
dev.off()




