library(inactiveXX)

devtools::install_github('constantAmateur/inactiveXX')


bams10X <- list.files('ASEReadCounter/individuals/sample', pattern='-?.bam')

# Step 1 - call heterozygous SNPS
hSNPs = hetSNPsFromRNA(bams10X,refGenome10X,outputs=sprintf('%s_scRNA_1kSNPs_XCnts.tsv',names(bams10X)))

# Step 2 - Get counts from heterozygous SNPs on X
XCnts = filterCountsX(hSNPs)

# Step 3 - Infer Xi status
fit = inferInactiveX(XCnts,nParallel=8)

#Visualise the fit
plotSolutitons(fit)

# Deviations in Xi skew within one individual
srat <- GetAssayData(pbmc)
milo = estCelltypeSkew(srat,fit,resultsPassthrough='annot')
plotDAbeeswarm(milo$res,group.by='annot')