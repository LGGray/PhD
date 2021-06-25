library(MAST)
library(Seurat)

# Read in files
pool_1 <- read.delim("~/datasets/OneK1k/M_vs_F_DEout/fcHurdleSig/pool_1.txt", sep = "\t", header = T)
pool_2 <- read.delim("~/datasets/OneK1k/M_vs_F_DEout/fcHurdleSig/pool_2.txt", sep = "\t", header = T)
pool_3 <- read.delim("~/datasets/OneK1k/M_vs_F_DEout/fcHurdleSig/pool_3.txt", sep = "\t", header = T)

# Read in all possible genes
genes <- read.delim("~/datasets/OneK1k/M_vs_F_DEout/genes.txt", sep="\t", header = F)[-1,]

# Creating function
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
colnames(outfile)[2] <- "combined_logFC"

# pull logFC and save to outfile
for (gene in genes){
  gene <- paste("^", gene,"$",sep = "")
  pval <- c(pool_1[grep(gene,pool_1[,1]),2], pool_2[grep(gene,pool_2[,1]),2], pool_3[grep(gene,pool_3[,1]),2])
  pval <- pval[!is.na(pval)]
  outfile[grep(gene, outfile),2] <- sumPvalsMethod(pval)
}

outfile <- outfile[!is.na(outfile$combined_logFC),]

write.table(outfile, file="~/datasets/OneK1k/M_vs_F_DEout/DEGs.txt", sep="\t", row.names = F)

