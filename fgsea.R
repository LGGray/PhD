library(fgsea)
library(msigdbr)
library(ggplot2)
library(data.table)

#msigdbr_df <- msigdbr(species = "human", category = "C7", subcategory = "IMMUNESIGDB")
msigdbr_df <- msigdbr(species = "human", category = "H")
pathways = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

setwd("~/datasets/OneK1k/C_vs_RA_DEout/nebula/")
files <- list.files(pattern = "\\.RDS$")
for (i in 1:length(files)){
  cell <- gsub(".RDS", "", files[i])
  print(cell)
  re <- readRDS(files[i])
  results.ord <- re$summary[ order(re$summary[,"logFC_ccY"]), ]
  ranks <- results.ord$p_ccY
  names(ranks) <- results.ord$gene
  ranks <- log2(ranks)
  
  fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize = 500)
  fwrite(fgseaRes, file=paste0("../fgsea/",cell,".H.tsv"), sep="\t", sep2=c("", " ", ""))
}


files = list.files("../fgsea/", pattern = "*.H.tsv")
plot_list = list()
for (i in 1:length(files)){
  cell <- gsub(".tsv", "", files[i])
  print(cell)
  fgseaRes <- read.delim(paste0("../fgsea/", files[i]))
  sigpath <- subset(fgseaRes, pval <= 0.05)
  sigpath.ord <- sigpath[order(sigpath[,"pval"]),]
  sigpath$pathway <- gsub("^GSE\\d+_", "", sigpath$pathway)
  
  p <- ggplot(sigpath[1:20,], aes(x = pval, y = pathway)) + 
    geom_point(aes(size = size, fill = pval), alpha = 0.75, shape = 21, stroke = 0.5) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(paste0(cell, "", " - Top 20 enriched H pathways")) +
    xlab("Pvalue")
  
  plot_list[[i]] = p
  
}

pdf("~/plots/fgsea/RA.H.cells.pdf", width = 14)
for (i in 1:length(files)){
  print(plot_list[[i]])
}
dev.off()
