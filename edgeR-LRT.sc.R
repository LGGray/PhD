library(edgeR)
library(Seurat)
library(reshape2)
library(ggplot2)
library(Matrix)
library(palettetown)
library(ggrepel)

args = commandArgs(trailingOnly=TRUE)
args = as.numeric(args[1])

pbmc <- readRDS("~/datasets/OneK1k/C_vs_RA_DEout/pbmc.subset.RDS")
# Select cell type
cell = levels(pbmc)[args]
print(cell)
# subset object by cell type
pbmc.cell <- subset(pbmc, predicted.celltype.l2 %in% cell)
expr <- GetAssayData(object = pbmc.cell, slot = "scale.data")
expr = as.data.frame(expr)
pbmc.cell$pool <- factor(pbmc.cell$pool)
pbmc.cell$RA <- factor(pbmc.cell$RA)
pbmc.cell$individual <- factor(pbmc.cell$individual)
# edgeR-LRT
targets = data.frame(group = pbmc.cell$RA,
                     pool = pbmc.cell$pool,
                     age = pbmc.cell$age,
                     individual = pbmc.cell$individual)
design <- model.matrix(~0+pool+age+group, data=targets)
y = DGEList(counts = expr, group = targets$group)
# Filter for expression in 5% of cells
keep <- rowSums(expr > 0) > dim(y)[2]*0.05
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design)
lrt <- glmLRT(fit, coef=dim(design)[2])
print(summary(decideTests(lrt)))
lrt <- as.data.frame(lrt)
lrt$FDR <- p.adjust(lrt$PValue, "BH")
gene <- rownames(lrt)
lrt <- cbind(gene,lrt)
write.table(lrt, paste0("~/datasets/OneK1k/C_vs_RA_DEout/edgeR-LRT/sc/", cell, ".txt"),
            row.names=F, sep="\t", quote = F)


setwd("~/datasets/OneK1k/C_vs_RA_DEout/edgeR-LRT/sc/")
# read in XCI escape data
load("~/datasets/OneK1k/X_escape/escapees.Rdata")
xcape_genes <- rownames(escape)

# Reading in edgeR-LRT output files
files <- list.files(pattern = "\\.txt$")

# Set ggrepel
options(ggrepel.max.overlaps = Inf)
# Create empty plot list
plot_list = list()
# Iterate over files and create volcano plot
for (i in 1:length(files)){
  cell <- gsub(".txt", "", files[i])
  print(cell)
  # Read in file
  res <- read.delim(files[i], header=T)
  if (dim(subset(res, FDR <= 0.05))[1] > 0){
    # Indicate significance
    threshold <- res$FDR <= 0.05
    res$threshold <- threshold
    # Order on p value
    res_ordered <- res[order(res$FDR),]
    # Indicate if gene is XCI escaping and DE
    res$xcape <- ifelse(res$gene %in% xcape_genes &
                          res$threshold == T, T, F)

    p = ggplot(res) +
      geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
      scale_color_poke(pokemon=155, spread=2) +
      geom_text_repel(box.padding = 0.75, aes(x=logFC, y=-log10(FDR), label = ifelse(threshold == T , gene, ""))) +
      ggtitle(paste(cell, "Volcano Plot")) +
      xlab("log2 fold change") +
      ylab("-log10 FDR") +
      theme(legend.position = "none",
            plot.title = element_text(size = rel(1.5), hjust = 0.5),
            axis.title = element_text(size = rel(1.25)))
    
    plot_list[[i]] = p

    # if (dim(subset(res, xcape == T))[1] > 0){
    #   print(T)
    #   # q = p + geom_point(aes(logFC[xcape==T], -log10(FDR)[xcape==T]), colour="black", shape=1, size=4)
    #   # plot_list[[i]] = q
    # } else {
    #   plot_list[[i]] = p
    # }
  }
}

pdf("~/plots/RA_genes/edgeR.sc.VP.pdf")
for (i in 1:length(files)){
  print(plot_list[[i]])
}
dev.off()