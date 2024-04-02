### Code to plot a figure of GWAS variants across sutoimmune diseases studied ###

library(reshape2)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)

setwd('/directflow/SCCGGroupShare/projects/lacgra/gwas/')

# Load data
CD <- read.delim('CD.tsv')
pSS <- read.delim('pSS.tsv')
SLE <- read.delim('SLE.tsv')
MS <- read.delim('MS.tsv')
UC <- read.delim('UC.tsv')

study_colours <- list(CD='#B00074', pSS='#B09F99', SLE='#6162B0', MS='#67B05D', UC='#187362')

variants <- lapply(list(CD, pSS, SLE, MS, UC), function(x) {
    tmp <- subset(x, pValue < 5e-8)$riskAllele
    # gsub('-.', '', tmp)
})
names(variants) <- c('CD', 'pSS', 'SLE', 'MS', 'UC')

pdf('gwas_upset.pdf', onefile=F)
upset(fromList(variants), order.by = 'freq', nsets=5, sets.bar.color = unlist(study_colours))
dev.off()

variants.chrX <- lapply(list(CD, pSS, SLE, MS, UC), function(x) {
    tmp <- subset(x, pValue < 5e-8)
    tmp[grep('X:', tmp$locations), 'riskAllele']
})
names(variants.chrX) <- c('CD', 'pSS', 'SLE', 'MS', 'UC')

pdf('gwas_upset_chrX.pdf', onefile=F)
upset(fromList(variants.chrX), order.by = 'freq', nsets=5)
dev.off()

variant.matrix <- fromList(variants)
rownames(variant.matrix) <- unique(unlist(variants))
sort(rowSums(variant.matrix))

pdf('gwas_heatmap.pdf', onefile=F)
Heatmap(as.matrix(variant.matrix), name = 'GWAS Variants', col = colorRamp2(c(0, 1), c('white', 'red')),
show_row_names = F, column_names_gp = gpar(fontsize = 12, angle = 0))
dev.off()
