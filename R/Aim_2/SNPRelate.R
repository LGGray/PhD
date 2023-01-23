# Dimension reduction of X chromosome genotypes

setwd("~/datasets/OneK1k")
library("SNPRelate")
library(ggplot2)
library(ggrepel)

# Convert vcf to gds file
vcf.fn <- '~/datasets/OneK1k/female.recode.vcf'
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- snpgdsOpen("ccm.gds")

# Remove SNPs in LD
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, autosome.only=F)
snpset.id <- unlist(unname(snpset))

ccm_pca <- snpgdsPCA(genofile, snp.id=snpset.id, autosome.only=F)

# variance proportion (%)
pc.percent <- ccm_pca$varprop*100
head(round(pc.percent, 2))

colnames(ccm_pca$eigenvect) <- paste0('PC', 1:32)

pdf('chrX.PCA.pdf')
ggplot(data.frame(ccm_pca$eigenvect), aes(x=PC2, y=PC1, label = ccm_pca$sample.id)) +
  geom_point() +
  geom_text_repel()
dev.off()

