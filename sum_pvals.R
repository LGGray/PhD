library(MAST)
library(Seurat)

# Read in files
setwd("~/datasets/OneK1k/M_vs_F_DEout/fcHurdleSig/")
fnames <- list.files()
df <- lapply(fnames, read.delim)

# Read in all possible genes
genes <- read.delim("~/datasets/OneK1k/M_vs_F_DEout/genes.txt", sep="\t", header = F)[-1,]

# Creating fishers method function
psumunif = function(x,n) 1/factorial(n) * sum(sapply(0:n, 
           function(k) (-1)^k * choose(n,k) * ifelse(x > k,x-k,0)^(n)))
sumPvalsMethod = function(x) {
  n = length(x)
  if (n < 10) {
    psumunif(sum(x),n)
  } else {
    pnorm(sum(x),n/2,sqrt(n/12),lower=TRUE)
  }
  }

# Create outfile dataframe
outfile <- as.data.frame(matrix(nrow=length(genes), ncol=2))
outfile[,1] <- genes
colnames(outfile)[1] <- "gene"
colnames(outfile)[2] <- "FDR"

# pull logFC and save to outfile
pval <- c()
for (gene in genes){
  gene <- paste0("^", gene,"$")
  for ( i in 1:length(df)){
    pval <- c(pval, df[[i]][grep(gene,df[[i]][,1]),2])
    if (i == length(df)){
      pval <- pval[!is.na(pval)]
      pval <- p.adjust(pval, method="fdr")
      outfile[grep(gene, outfile[,1]),2] <- sumPvalsMethod(pval)
    }
  }
}


outfile <- subset(outfile, FDR <= 0.01)

write.table(outfile, file="~/datasets/OneK1k/M_vs_F_DEout/DEGs.txt", sep="\t", row.names = F)
