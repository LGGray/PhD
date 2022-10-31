library(clusterProfiler)
source('~/R_code/functions/edgeR.list.R')
source('~/R_code/functions/chisq.test.degs.R')

gene.set <- read.gmt('~/gene.sets/msigdb.v7.5.1.symbols.gmt')

gene.set <- data.frame(gene.set)

ER <- gene.set[grep('ESTROGEN', gene.set$term),]
AR <- gene.set[grep('ANDROGEN', gene.set$term),]
PR <- gene.set[grep('PROGESTERONE', gene.set$term),]

edgeR <- edgeR.list('/directflow/SCCGGroupShare/projects/lacgra/sex.bias/edgeR-LRT/', filter=F)



ER.enrich <- lapply(edgeR, function(x){chisq.test.degs(x, unique(ER$gene), logFC=0.5)$p.value})
AR.enrich <- lapply(edgeR, function(x){chisq.test.degs(x, unique(AR$gene), logFC=0.5)$p.value})
PR.enrich <- lapply(edgeR, function(x){chisq.test.degs(x, unique(PR$gene), logFC=0.5)$p.value})

ER.enrich[sapply(ER.enrich, function(x) x < 0.05)]
AR.enrich[sapply(AR.enrich, function(x) x < 0.05)]
PR.enrich[sapply(PR.enrich, function(x) x < 0.05)]


ER.DC <- subset(edgeR$DC, FDR < 0.05 & abs(logFC) > 0.5 & gene %in% unique(ER$gene))
PR.DC <- subset(edgeR$DC, FDR < 0.05 & abs(logFC) > 0.5 & gene %in% unique(PR$gene))
PR.NK_CD56bright <- subset(edgeR$NK_CD56bright, FDR < 0.05 & abs(logFC) > 0.5 & gene %in% unique(PR$gene))

setwd('/directflow/SCCGGroupShare/projects/lacgra/sex.bias/')
write.table(ER.DC, 'ER.DC.txt', row.names=F, quote=F)
write.table(PR.DC, 'PR.DC.txt', row.names=F, quote=F)
write.table(PR.NK_CD56bright, 'PR.NK_CD56bright.txt', row.names=F, quote=F)

#-----------------
female <- read.delim('sex_specific_female_eQTLs_not_in_joint_results_and_passed_both_tests.tsv')

female <- split(female, female$celltype)
ER.ora <- lapply(female, function(x){
  enricher(gene = x$geneid,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.01,
           TERM2GENE = ER)
})

PR.ora <- lapply(female, function(x){
  enricher(gene = x$geneid,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.01,
           TERM2GENE = PR)
})

AR.ora <- lapply(female, function(x){
  enricher(gene = x$geneid,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.01,
           TERM2GENE = AR)
})

male <- read.delim('sex_specific_male_eQTLs_not_in_joint_results_and_passed_both_tests.tsv')
male <- split(male, male$celltype)
ER.ora <- lapply(male, function(x){
  enricher(gene = x$geneid,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.01,
           TERM2GENE = ER)
})

PR.ora <- lapply(male, function(x){
  enricher(gene = x$geneid,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.01,
           TERM2GENE = PR)
})

AR.ora <- lapply(male, function(x){
  enricher(gene = x$geneid,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.01,
           TERM2GENE = AR)
})

lapply(AR.ora, nrow)
ER.ora[lapply(PR.ora, function(x) nrow(x) > 0)]

data.frame(AR.ora$CD8_TEM)
