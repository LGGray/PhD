library(dplyr)
library(ggplot2)

dir <- '/directflow/SCCGGroupShare/projects/lacgra/'
setwd(paste0(dir,'ASEReadCounter/individuals/'))

individuals <- read.delim(paste0(dir,'/datasets/OneK1k/female.ids.txt'), header=F)$V1
individuals <- c('684_685', '685_686', '687_688', '692_693', '693_694')

load("/directflow/SCCGGroupShare/projects/sarba2/data/onek1k/samples/metadata.Rdata")
metadata$cell.id <- gsub('-\\d', '', rownames(metadata))
keep <- unique(metadata$predicted.celltype.l3)[-c(9,29)]
metadata <- subset(metadata, predicted.celltype.l3 %in% keep)

# Read in files
ASE <- lapply(individuals[1:5], function(x){
  files <- list.files(x, pattern='-1.txt', full.names = T)
  tmp <- lapply(files, function(x){
    read.delim(x, colClasses= c(rep('character',5), rep('numeric', 8)))
  })
  names(tmp) <- as.character(gsub('-1.txt', '', basename(files)))
  tmp[sapply(tmp, function(x) nrow(x) > 1)]
})
names(ASE) <- individuals[1:5]

# Bind all cells into a df
ASE.df <- lapply(ASE, function(x){
  dplyr::bind_rows(x, .id='cell.id') %>%
  left_join(., metadata[,c(23, 16)], by='cell.id')
})

# Sum counts for each SNP within each celltype
ASE.sum <- lapply(ASE.df, function(x){
  x %>% 
    group_by(position, predicted.celltype.l3) %>% 
    summarise(refCount = sum(refCount), altCount = sum(altCount), totalCount = sum(totalCount)) %>%
    mutate(position=as.numeric(position)) %>%
    mutate(hap1 = refCount/totalCount, hap2 = altCount/totalCount) %>%
    mutate(max = round(pmax(hap1, hap2), 2)) %>%
    mutate(coord=paste('X', position, position+1, sep=':')) %>%
    filter(predicted.celltype.l3 != 'NA') %>%
    as.data.frame()
})
unique(ASE.sum[[1]]$predicted.celltype.l3)
save(ASE.sum, file='analysis/ASE.sum.Rdata')

# Filter summed counts by nonPAR
ASE.sum.nonPAR <- lapply(ASE.sum, function(x){
  x[x$position > 2699520 & x$position < 154931044,]
})
save(ASE.sum.nonPAR, file='ASE.sum.nonPAR.Rdata')


ggplot(ASE.sum.nonPAR[[3]], aes(x = max)) +
  facet_wrap(~predicted.celltype.l3) +
  geom_histogram()
fit <- MASS::fitdistr(ASE.sum.nonPAR[[3]]$max, "beta-binomial")
print(fit)

# load the fitdistrplus package and VGAM
library(fitdistrplus)

# fit a beta-binomial distribution to the data
fit <- fitdistr(ASE.sum.nonPAR[[1]]$max, distr="normal", start=c(0.5, 100))

# print a summary of the fit
data <- ASE.sum.nonPAR[[1]]$max
fit <- fitdist(data, distr = "betabinom", start=list(prob=0.5, size=1875), method="mge")
summary(fit)

# plot the fit
plot(fit)


# Heatmap
ASE.sum[[1]]$scaled <- scale(ASE.sum[[1]]$max)
mtx <- dcast(ASE.sum[[1]], position ~ predicted.celltype.l3, value.var = 'scaled') %>%
  tibble::column_to_rownames('position') %>%
  replace(is.na(.), 0)
m <- c( "average", "single", "complete", "ward")
ac <- function(x) {
  agnes(mtx, method = x)$ac
}
clust <- agnes(t(mtx), method = m[which.max(unlist(lapply(m, ac)))])
pltree(clust, cex = 0.6, hang = -1, main = "Dendrogram of agnes")

ggplot(ASE.sum[[1]], aes(x=factor(predicted.celltype.l3, level=colnames(mtx)[clust$order]),
                         y=as.character(position), fill=scaled)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="purple", high="red") +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.text.y = element_blank())

pdf('analysis/test.nonPAR.distribution.pdf')
tmp <- lapply(df, function(x){
  hist(x$max, main = unique(x$predicted.celltype.l3))
})
dev.off()


ASE.sum.nonPAR
ggplot(ASE.sum.nonPAR[[1]], aes(x=max, fill=predicted.celltype.l3)) +
  geom_histogram(position='stack') +
  facet_grid(predicted.celltype.l3 ~ ., scales='free')



qqnorm(ASE.sum.nonPAR[[3]]$max)
qqline(ASE.sum.nonPAR[[3]]$max)

x <- ASE.sum.nonPAR[[3]]$max
y <- gamlss.dist::rBB(n=length(x), mu=0.5)
ks.test(x,y)
plot(ecdf(x=x))
lines(ecdf(x = y), col = 2)

b <- gamlss.dist::rBB(n=length(x), mu=0.5)
plot(density(x))

  

library(biomaRt)

listEnsemblArchives()

snp_mart <- useEnsembl(biomart="ENSEMBL_MART_SNP", 
                       host="grch37.ensembl.org", 
                       dataset="hsapiens_snp")

snp.id <- getBM(attributes = c('refsnp_id','chrom_start', 'associated_gene'), 
      filters = 'chromosomal_region', 
      values = ASE.sum.nonPAR[[1]]$coord[1:10], 
      mart = snp_mart)


