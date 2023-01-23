library(EGAD)
library(Seurat)
library(Biobase)
load('~/tools/CoExpNets/bin/run_GBA.Rdata')
source('~/tools/CoExpNets/bin/helper_functions.r')
load('~/tools/CoExpNets/data/GO.human.Rdata')
data(biogrid)
library(ggplot2)

# Read in file
pbmc <- readRDS("datasets/OneK1k/C_vs_RA_DEout/pbmc.RDS")

# X chromosome genes nonPAR
load("~/datasets/OneK1k/X_escape/chrX.Rdata")
#load("datasets/OneK1k/X_escape/escapees.Rdata")

# DEGs for celltype
deg <- readRDS("~/datasets/OneK1k/C_vs_RA_DEout/nebula/Treg.RDS")
degX <- subset(deg$summary, p_ccY <= 0.05 & gene %in% rownames(chrX))

# Set empty matrix
nrows <- dim(degX)[1]
aggs = diag(nrows)

# select individual ids 
ids <- unique(pbmc$individual)

# open plot list
plot_list = list()

# Construct network for each sample 
for (i in 1:length(ids)){
  # subset for an individual
  data <- subset(pbmc, subset = individual %in% ids[i] & 
                   predicted.celltype.l2 == "Treg")
  RA <- unique(data$RA)
  
  # extract scaled data from subset
  df <- t(FetchData(data, vars = degX$gene))
  
  # contruct network
  set <- new("ExpressionSet", exprs=as.matrix(df))
  net <- build_coexp_expressionSet(set, rownames(df))
  med <- median(net, na.rm = T)
  net[is.na(net)] = med
  
  save(net, file=paste0("~/datasets/OneK1k/C_vs_RA_DEout/co-expression/Treg.", ids[i], ".", RA, ".Rdata"))
  
  net <- melt(net)
  net$Var1 <- gsub("rna_", "", net$Var1)
  net$Var2 <- gsub("rna_", "", net$Var2)
  
  # construct heatmap
  p <- ggplot(data = net, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + 
    scale_fill_gradient2(low = "white", high = "red", mid = "blue", 
                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                         name = "Correlation") +
    ggtitle(paste("Treg", ids[i], RA, "Co-expression matrix", sep = " ")) +
    xlab("Gene") + ylab("Gene") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1)) +
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = round(value, 2)), color = "black", size = 4)
  
  plot_list[[i]] = p
}

pdf("~/plots/co-express/Treg.pdf")
for (i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()

setwd("~/datasets/OneK1k/C_vs_RA_DEout/co-expression")
files <- list.files(pattern = paste0("Y.Rdata"))
  for ( i in 1:length(files)){
    load(files[i])
    if (i == 1) {
      agg = net
    } else {
      agg = net + agg
    }
  }
  
agg.rank = matrix(rank(agg, na.last = "keep", ties.method = "average"), 
              nrow=dim(agg)[1], ncol=dim(agg)[2] )
rownames(agg.rank) = rownames(agg)
colnames(agg.rank) = rownames(agg)
agg.rank = agg.rank/max(agg.rank, na.rm=T)

save(agg.rank, file = paste0("Treg.", "Y", ".Rdata"))

agg.rank <- melt(agg.rank)
pdf(paste0("~/plots/co-express/Treg.", "Y", ".agg.pdf"))
ggplot(data = agg.rank, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "white", high = "red", mid = "blue", 
                     midpoint = 0.5, limit = c(0,1), space = "Lab", 
                     name="Correlation") +
  ggtitle(paste("Treg", "Y", "Aggregated Co-expression matrix", sep=" ")) +
  xlab("Gene") + ylab("Gene") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                 size = 12, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = round(value, 2)), color = "black", size = 4)
dev.off()


