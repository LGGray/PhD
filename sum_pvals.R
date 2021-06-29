library(MAST)
library(Seurat)

# Read in files
setwd("~/datasets/OneK1k/M_vs_F_DEout/fcHurdleSig/")
fnames <- list.files()
df <- lapply(fnames, read.delim)

# Read in all possible genes
genes <- read.delim("~/datasets/OneK1k/M_vs_F_DEout/genes.txt", sep="\t", header = F)[-1,]

# Methods of combining p-values:  
fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)
psumunif = function(x,n) 1/factorial(n) * sum(sapply(0:n, function(k) (-1)^k * choose(n,k) * ifelse(x > k,x-k,0)^(n)))
sumPvalsMethod = function(x) {
  n = length(x)
  if (n < 10) {
    psumunif(sum(x),n)
  } else {
    pnorm(sum(x),n/2,sqrt(n/12),lower=TRUE)
  }
}
binomMethod = function(x) pbinom(sum(x < .05),size=length(x),prob=.05,lower=FALSE)

# Create outfile dataframe
outfile <- as.data.frame(matrix(nrow=length(genes), ncol=4))
outfile[,1] <- genes
colnames(outfile)[1] <- "genes"
colnames(outfile)[2] <- "Pval - sum"
colnames(outfile)[3] <- "Pval - binom"
colnames(outfile)[4] <- "Pval - fishers"

# pull pvals, adjust then save to outfile
pval <- c()
for (gene in genes){
  gene <- paste0("^", gene,"$")
  for ( i in 1:length(df)){
    pval <- c(pval, df[[i]][grep(gene,df[[i]][,1]),2])
    if (i == length(df)){
      pval <- pval[!is.na(pval)]
      # Need to do this across all the p-values at once, or you need to specify the number of p-values you are looking at  
      # pval <- p.adjust(pval, method="fdr")
      outfile[grep(gene, outfile[,1]),2] <- sumPvalsMethod(pval)
      outfile[grep(gene, outfile[,1]),3] <- binomMethod(pval)
      outfile[grep(gene, outfile[,1]),4] <- fishersMethod(pval)
    }
  }
}

# Do the FDR calc here 
padjusted = cbind(p.adjust(outfile[,2]), p.adjust(outfile[,3]), p.adjust(outfile[,4]))  
outfile = cbind(outfile, padjusted)

write.table(outfile, file="~/datasets/OneK1k/M_vs_F_DEout/DEGs.txt", sep="\t", row.names = F)
