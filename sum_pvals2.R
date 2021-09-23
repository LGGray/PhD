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

# read in RDS file containing MAST output
obj = readRDS("/directflow/SCCGGroupShare/projects/lacgra/fcHurdleSig_object.RDS")

# order each object by primerid alphabetically 
o1 = sapply(1:length(obj), function(i) order(obj[[i]][,1]) ) 

# pull out the pvalue and fdr from each file
pvals = sapply(1:length(obj), function(i) obj[[i]]$pvalue[o1[,i]])
fdrss = sapply(1:length(obj), function(i) obj[[i]]$fdr[o1[,i]]) 

# perform function
pvals1 = sapply(1:dim(pvals)[1], function(i) sumPvalsMethod(pvals[i,]))
pvals2 = sapply(1:dim(pvals)[1], function(i) binomMethod(pvals[i,]))
pvals3 = sapply(1:dim(pvals)[1], function(i) fishersMethod(pvals[i,]))
pvals1f = sapply(1:dim(pvals)[1], function(i) sumPvalsMethod(fdrss[i,]))
pvals2f = sapply(1:dim(pvals)[1], function(i) binomMethod(fdrss[i,]))
pvals3f = sapply(1:dim(pvals)[1], function(i) fishersMethod(fdrss[i,]))

# Create empty dataframe and fill with data
outfile <- as.data.frame(matrix(nrow=length(pvals1), ncol=7))
outfile[,1] <- obj[[1]][o1[,1],1]
outfile[,2] <- p.adjust(pvals1)
outfile[,3] <- p.adjust(pvals2)
outfile[,4] <- p.adjust(pvals3)
outfile[,5] <- p.adjust(pvals1f)
outfile[,6] <- p.adjust(pvals2f)
outfile[,7] <- p.adjust(pvals3f)

# Replace column names
colnames(outfile) <- c("primerid", "pval-sum", "pval-binom", "pval-fishers", "fdr-sum", "fdr-binom", "fdr-fishers")

# Save output
write.table(outfile, "datasets/OneK1k/M_vs_F_DEout/DEGs.txt", sep="\t", row.names=F)
