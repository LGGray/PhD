library(EGAD)
library(Seurat)
library(Biobase)
load('~/tools/CoExpNets/bin/run_GBA.Rdata')
source('~/tools/CoExpNets/bin/helper_functions.r')
load('~/tools/CoExpNets/data/GO.human.Rdata')
data(biogrid)

args = commandArgs(trailingOnly=TRUE)

# Create directories to store results
mainDir <- '~/datasets/OneK1k/co-express/h_net'
subDir <- paste0('pool_', args[1])
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

mainDir <- '~/datasets/OneK1k/co-express/networks'
subDir <- paste0('pool_', args[1])
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

mainDir <- '~/plots/co-express'
subDir <- paste0('pool_', args[1])
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

# read in data
pbmc <- readRDS(paste0('~/datasets/OneK1k/M_vs_F_DEout/pbmc_RDS/pbmc', args[1],'.RDS'))
print(pbmc)
degs <- read.delim("~/datasets/OneK1k/M_vs_F_DEout/DEGs.txt", sep="\t", header = T)

# read in escape data
load("/home/lacgra/datasets/OneK1k/X_escape/de_sex_genes.Rdata")
chrx <- subset(annotations, annotations[,2] == 1 & annotations[,3] == 0)
escape <- subset(annotations, annotations[,19] == 1 | annotations[,20] == 1)
chrX <- subset(degs, gene %in% rownames(chrx))
escapees <- subset(degs, gene %in% rownames(escape))

# select individual ids 
ids <- unique(pbmc$individual)
ids_sex <- unique(paste0(pbmc$individual, "_", pbmc$sex))

# Set empty matrix
nrows <- dim(GetAssayData(object = pbmc)[1])
aggs = diag(nrows)

# Aggregate 
for (i in 1:length(ids)){
  # subset for an individual
  data <- subset(pbmc, subset = individual %in% ids[i])
  
  # extract scaled data from subset
  df <- t(FetchData(data, vars = chrX[,1]))
  
  # contruct network
  set <- new("ExpressionSet", exprs=as.matrix(df))
  net <- build_coexp_expressionSet(set, rownames(df))
  med <- median(net, na.rm = T)
  net[is.na(net)] = med
  
  save(net, file=paste0("~/datasets/OneK1k/co-express/networks/", 
                        paste0('pool_', args[1], '/'),"net_", i, ".", "Rdata"))
  
  # construct heatmap
  pdf(paste0("~/plots/co-express/", paste0('pool_', args[1], '/'), "network_heatmap.", i, ".pdf"))
  x <- heatmap.2(net, density="none", trace="none")
  dev.off()
  save(x, file=paste0("~/datasets/OneK1k/co-express/h_net/", paste0('pool_', args[1], '/'), 'h_net.', i, ".Rdata"))
  print(i)
  
  if (i == 1) {
    agg = net
  } else {
    agg = net + agg
  }
}

agg.rank =  matrix(rank(agg, na.last = "keep", ties.method = "average"), nrow=dim(agg)[1], ncol=dim(agg)[2] )
rownames(agg.rank) = rownames(agg)
colnames(agg.rank) = rownames(agg)
agg.rank = agg.rank/max(agg.rank, na.rm=T)

save(agg.rank, file = paste0("~/datasets/OneK1k/co-express/networks/", paste0('pool_', args[1], '/'), "coexp.agg.rank.Rdata"))

# set annotations
genelist <- rownames(agg.rank)
goterms <- unique(GO.human[,3])
annotations <- make_annotations(GO.human[,c(1,3)],genelist,goterms)

# Assess output
aurocs = run_GBA(agg.rank, annotations)
aurocs[[2]] = ""
save(agg.rank, file = paste0("~/datasets/OneK1k/co-express/networks/", paste0('pool_', args[1], '/'), "aurocs.Rdata"))
pdf(paste0("~/plots/co-express/", paste0('pool_', args[1], '/'), "aurocs.pdf"))
plot_density_compare(aurocs[[1]][,1], aurocs[[1]][,3])
dev.off()
