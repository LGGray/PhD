# Code to analyse differentially expressed genes in the ChEA3 database

# Read in DEG files
files <- list.files('datasets/integrated/psuedobulk/downregulated/', full.names = T)
down <- lapply(files, function(x){
  df <- read.delim(x)
})
names(down) <- gsub('.txt', '', basename(files))

# Access ChEA3 database from R
library(httr)
library(jsonlite)

genes = down[[1]]$gene

url = "https://maayanlab.cloud/chea3/api/enrich/"
encode = "json"
payload = list(query_name = "myQuery", gene_set = genes)

#POST to ChEA3 server
response = POST(url = url, body = payload, encode = encode)
json = content(response, "text")

#Integrated--meanRank results as list of R dataframes
results = fromJSON(json)[[1]]

# Load gene sets
load('datasets/XCI/chrX.Rdata')
load('datasets/XCI/escapees.Rdata')
immport <- read.delim('gene.sets/Immport.GeneList.txt')

# Create list of genes targeted by TF
TF <- lapply(results$Overlapping_Genes, function(x){
  unlist(strsplit(x, ','))
})
names(TF) <- results$TF

immune.TF <- TF[names(TF) %in% immport$Symbol]

# Which TFs escape XCI?
xcape.TF <- TF[names(TF) %in% rownames(escape)]

# Which TFs target xcape genes?
xcape.target <- lapply(TF, function(x){
  x[x %in% rownames(escape)]
})
xcape.target <- xcape[sapply(xcape, function(x) length(x) > 0)]

# xcape.target by immune TF
xcape.target.immune <- xcape.target[names(xcape.target) %in% names(immune.TF)]

# DEGs regulated by oestradiol
lapply(TF[grep('ESR[1-2]', names(TF))], function(x){
  x[x %in% rownames(chrX)]
})
# DEGs regulated by androgen
lapply(TF[grep('^AR\\b', names(TF))], function(x){
  x[x %in% rownames(chrX)]
})
# DEGs regulated by progesterone
lapply(TF[grep('^PGR\\b', names(TF))], function(x){
  x[x %in% rownames(chrX)]
})
