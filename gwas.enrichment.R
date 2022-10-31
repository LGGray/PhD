source('~/external/ClusterHome/R_code/functions/edgeR.list.R')
source('~/external/ClusterHome/R_code/functions/chisq.test.degs.R')
load('~/external/ClusterHome/datasets/XCI/chrX.Rdata')

setwd('~/external/ClusterHome/datasets/integrated/psuedobulk/')

lof <- read.delim('../../../gwas/coding.variants.txt', header = F)
variants <- read.delim('../../../gwas/autoimmune.variants.txt', header = F)

all <- edgeR.list('.', filter = F)
deg <- edgeR.list('.', logfc = 0.5)
deg.lof <- lapply(deg, function(x){
  subset(x, gene %in% lof$V1)
})
deg.lof <- deg.lof[sapply(deg.lof, function(x) nrow(x) > 0)]
chi.lof <- lapply(all, function(x){
  chisq.test.degs(x, genes = lof$V1, logfc = 0.5)
})


deg.var <- lapply(deg, function(x){
  subset(x, gene %in% variants$V1)
})
deg.count <- lapply(deg.var, function(x){
  up <- subset(x, logFC > 0)
  down <- subset(x, logFC < 0)
  df.up <- data.frame(var=length(up$gene[up$gene %in% variants$V1]),total=length(up$gene), direction='Up')
  df.down <- data.frame(var=length(down$gene[down$gene %in% variants$V1]), total=length(down$gene), direction='Down')
  df <- rbind(df.up, df.down)
  df$total <- df$total-df$var
  return(df)
})
deg.count <- dplyr::bind_rows(deg.count, .id='celltype')
ggplot(melt(deg.count), aes(x=celltype, y=value, fill=direction)) + 
  geom_bar(stat='identity', position = position_dodge()) + 
  scale_fill_brewer(palette='Paired') +
  theme(axis.text.y = element_text(size=14)) +
  coord_flip() +
  xlab('') + ylab('DEG count') + labs(fill='Direction')


deg.var <- deg.var[sapply(deg.var, function(x) nrow(x) > 0)]
lapply(deg.var, nrow)
chi.var <- lapply(all, function(x){
  chisq.test.degs(x, genes = variants$V1, logfc = 0.5)
})

chi.stat <- lapply(chi.var, function(x) x$statistic)
pvalue <- lapply(chi.var, function(x) x$p.value)

df <- data.frame(stat=unlist(chi.stat), pvalue=unlist(pvalue))
df$cell <- gsub('.X-squared', '', rownames(df))
df$p.adjust <- p.adjust(df$pvalue)

ggplot(df, aes(x=stat, y=-log10(pvalue))) +
  geom_point() +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")

ggplot(df, aes(x=stat, y=-log10(p.adjust), label=cell)) +
  geom_point() +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
  geom_text(aes(label=ifelse(p.adjust<0.05,as.character(cell),'')), hjust=0, vjust=0,
            position=position_jitter(width = 1, height = 0.25))


