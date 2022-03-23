# Rscript to perform GSEA with fgsea
# msigdbr is used to download the gene pathways from MSigDB
# http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

library(fgsea)
library(msigdbr)
library(ggplot2)
library(data.table)

### MSigDB Hallmark gene sets
# msigdbr_df <- msigdbr(species = "human", category = "H")
### MSigDB immune gene sets
# msigdbr_df <- msigdbr(species = "human", category = "C7", subcategory = "IMMUNESIGDB")

pathways = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

# Set to directory with edgeR output. Example filename 'B_memory.txt'
setwd("")
files <- list.files(pattern = "\\.txt$")
for (i in 1:length(files)){
  cell <- gsub(".txt", "", files[i])
  print(cell)
  res <- read.delim(files[i], header=T)
  ordered <- res[order(res$logFC),]
  ranks <- ordered$logFC
  names(ranks) <- ordered$gene
  
  fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize = 500)
  fwrite(fgseaRes, file=paste0("output/path/",cell,".tsv"), sep="\t", sep2=c("", " ", ""))
}

# files = list.files("input/path/", pattern = ".tsv")
# plot_list = list()
# for (i in 1:length(files)){
#   cell <- gsub(".tsv", "", files[i])
#   print(cell)
#   fgseaRes <- read.delim(paste0("datasets/OneK1k/C_vs_RA_DEout/fgsea/all/", files[i]))
#   sigpath <- subset(fgseaRes, padj <= 0.05)
#   sigpath.ord <- sigpath[order(sigpath[,"padj"]),]
#   sigpath$pathway <- gsub("^GSE\\d+_", "", sigpath$pathway)
#   
#   p <- ggplot(sigpath[1:20,], aes(x = padj, y = pathway)) + 
#     geom_point(aes(size = size, fill = padj), alpha = 0.75, shape = 21, stroke = 0.5) +
#     scale_fill_gradient(low = "#e80200", high = "#284060") +
#     ggtitle(paste0(cell, "", " - Top 20 enriched IMMUNESIGDB pathways")) +
#     xlab("padj")
#   
#   plot_list[[i]] = p
#   
# }
# 
# pdf("~/plots/fgsea/all.C7.pdf", width = 14)
# for (i in 1:length(files)){
#   print(plot_list[[i]])
# }
# dev.off()
