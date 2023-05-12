
edgeR.files <- list.files('differential.expression/edgeR', full.names = T)
MAST.files <- list.files('differential.expression/MAST', full.names = T)

all.equal(basename(edgeR.files), basename(MAST.files))

edgeR <- lapply(edgeR.files, function(x){read.delim(x)})
names(edgeR) <- gsub('.txt', '', basename(edgeR.files))
MAST <- lapply(MAST.files, function(x){read.delim(x)})
names(MAST) <- gsub('.txt', '', basename(MAST.files))

# Compare edgeR and MAST
intersection <- lapply(1:length(edgeR), function(x){
    merge <- merge(edgeR[[x]], MAST[[x]], by='gene')
    merge.deg <- subset(merge, FDR < 0.05 | p_val_adj < 0.05)
    merge.deg
})
intersection[[2]][,1:3]

# Calculate enrichment of chrX genes
source('../../PhD/functions/chisq.test.degs.R')
load('../../datasets/XCI/chrX.Rdata')

edgeR.enrichment <- lapply(edgeR, function(x){chisq.test.edgeR(x, rownames(chrX), 0.1)})
names(edgeR.enrichment) <- names(edgeR)
MAST.enrichment <- lapply(MAST, function(x){chisq.test.MAST(x, rownames(chrX), 0.1)})
names(MAST.enrichment) <- names(MAST)
