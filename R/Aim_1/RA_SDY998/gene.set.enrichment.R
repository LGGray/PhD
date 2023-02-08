library(ggplot2)
library(reshape2)
library(UpSetR)
library(clusterProfiler)
library(dplyr)

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')

# Load Gene Sets
kegg <- read.gmt('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/c2.cp.kegg.v7.5.1.symbols.gmt')
reactome <- read.gmt('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/c2.cp.reactome.v7.5.1.symbols.gmt')
BP <- read.gmt('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/c5.go.bp.v7.5.1.symbols.gmt')
hallmark <- read.gmt('/directflow/SCCGGroupShare/projects/lacgra/gene.sets/h.all.v7.5.1.symbols.gmt')

setwd('/directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/RA_SDY998/')

# Load the unfiltered data
deg.list <- edgeR.list('psuedobulk', filter=F)
names(deg.list) <- gsub('.edgeR-LRT', '', names(deg.list))

# Over-representation analysis for downregulated genes: kegg
ora.down.kegg <- lapply(deg.list, function(x){
    result <- enricher(gene = subset(x, FDR < 0.05 & logFC < -0)$gene, 
           universe = subset(x, logFC < -0)$gene,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.05,
           TERM2GENE = kegg)
    if(!is.null(result) && nrow(data.frame(result)) > 0){
    return(subset(result@result, qvalue < 0.05))
  }
})
ora.down.kegg <- ora.down.kegg[!sapply(ora.down.kegg, is.null)]

# Over-representation analysis for upregulated genes: kegg
ora.up.kegg <- lapply(deg.list, function(x){
    result <- enricher(gene = subset(x, FDR < 0.05 & logFC > 0)$gene, 
           universe = subset(x, logFC > 0)$gene,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.05,
           TERM2GENE = kegg)
    if(!is.null(result) && nrow(data.frame(result)) > 0){
    return(subset(result@result, qvalue < 0.05))
  }
})
ora.up.kegg <- ora.up.kegg[!sapply(ora.up.kegg, is.null)]

# Over-representation analysis for downregulated genes: reactome
ora.down.reactome <- lapply(deg.list, function(x){
    result <- enricher(gene = subset(x, FDR < 0.05 & logFC < -0)$gene, 
           universe = subset(x, logFC < -0)$gene,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.05,
           TERM2GENE = reactome)
    if(!is.null(result) && nrow(data.frame(result)) > 0){
    return(subset(result@result, qvalue < 0.05))
  }
})
ora.down.reactome <- ora.down.reactome[!sapply(ora.down.reactome, is.null)]

# Over-representation analysis for upregulated genes: reactome
ora.up.reactome <- lapply(deg.list, function(x){
    result <- enricher(gene = subset(x, FDR < 0.05 & logFC > 0)$gene, 
           universe = subset(x, logFC > 0)$gene,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.05,
           TERM2GENE = reactome)
    if(!is.null(result) && nrow(data.frame(result)) > 0){
    return(subset(result@result, qvalue < 0.05))
  }
})
ora.up.reactome <- ora.up.reactome[!sapply(ora.up.reactome, is.null)]

# Over-representation analysis for downregulated genes: BP
ora.down.BP <- lapply(deg.list, function(x){
    result <- enricher(gene = subset(x, FDR < 0.05 & logFC < -0)$gene, 
           universe = subset(x, logFC < -0)$gene,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.05,
           TERM2GENE = BP)
    if(!is.null(result) && nrow(data.frame(result)) > 0){
    return(subset(result@result, qvalue < 0.05))
  }
})
ora.down.BP <- ora.down.BP[!sapply(ora.down.BP, is.null)]

# Over-representation analysis for upregulated genes: BP
ora.up.BP <- lapply(deg.list, function(x){
    result <- enricher(gene = subset(x, FDR < 0.05 & logFC > 0)$gene, 
           universe = subset(x, logFC > 0)$gene,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.05,
           TERM2GENE = BP)
    if(!is.null(result) && nrow(data.frame(result)) > 0){
    return(subset(result@result, qvalue < 0.05))
  }
})
ora.up.BP <- ora.up.BP[!sapply(ora.up.BP, is.null)]

# Over-representation analysis for downregulated genes: hallmark
ora.down.hallmark <- lapply(deg.list, function(x){
    result <- enricher(gene = subset(x, FDR < 0.05 & logFC < -0)$gene, 
           universe = subset(x, logFC < -0)$gene,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.05,
           TERM2GENE = hallmark)
    if(!is.null(result) && nrow(data.frame(result)) > 0){
    return(subset(result@result, qvalue < 0.05))
  }
})
ora.down.hallmark <- ora.down.hallmark[!sapply(ora.down.hallmark, is.null)]

# Over-representation analysis for upregulated genes: hallmark
ora.up.hallmark <- lapply(deg.list, function(x){
    result <- enricher(gene = subset(x, FDR < 0.05 & logFC > 0)$gene, 
           universe = subset(x, logFC > 0)$gene,
           pAdjustMethod = "fdr", 
           qvalueCutoff = 0.05,
           TERM2GENE = hallmark)
    if(!is.null(result) && nrow(data.frame(result)) > 0){
    return(subset(result@result, qvalue < 0.05))
  }
})
ora.up.hallmark <- ora.up.hallmark[!sapply(ora.up.hallmark, is.null)]

# Combine results for each cell type
ora.down.kegg <- bind_rows(ora.down.kegg, .id="celltype")
ora.up.kegg <- bind_rows(ora.up.kegg, .id="celltype")
ora.down.reactome <- bind_rows(ora.down.reactome, .id="celltype")
ora.up.reactome <- bind_rows(ora.up.reactome, .id="celltype")
ora.down.BP <- bind_rows(ora.down.BP, .id="celltype")
ora.up.BP <- bind_rows(ora.up.BP, .id="celltype")
ora.down.hallmark <- bind_rows(ora.down.hallmark, .id="celltype")
ora.up.hallmark <- bind_rows(ora.up.hallmark, .id="celltype")

# Create directory to save results
if(!dir.exists('ORA')){dir.create('ORA')}

# Save results
save(ora.down.kegg, file = 'ORA/ora.down.kegg.Rdata')
save(ora.up.kegg, file = 'ORA/ora.up.kegg.Rdata')
save(ora.down.reactome, file = 'ORA/ora.down.reactome.Rdata')
save(ora.up.reactome, file = 'ORA/ora.up.reactome.Rdata')
save(ora.down.BP, file = 'ORA/ora.down.BP.Rdata')
save(ora.up.BP, file = 'ORA/ora.up.BP.Rdata')
save(ora.down.hallmark, file = 'ORA/ora.down.hallmark.Rdata')
save(ora.up.hallmark, file = 'ORA/ora.up.hallmark.Rdata')

# function to split geneID column by '/' and find chrX genes in each
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/chrX.Rdata')
match.chrX <- function(x){
    if(!is.null(x)){
            result <- strsplit(x$geneID, split = "/")
        contains.chrX <- unlist(lapply(result, function(x) x %in% rownames(chrX)))
        if(any(contains.chrX)){
            out <- x[contains.chrX, ]
            return (out[!is.na(out$celltype),])
        } else {
            return (NULL)
        } 
    }
}

# Find chrX genes in each cell type
ora.down.kegg.chrX <- match.chrX(ora.down.kegg)
ora.up.kegg.chrX <- match.chrX(ora.up.kegg)
ora.down.reactome.chrX <- match.chrX(ora.down.reactome)
ora.up.reactome.chrX <- match.chrX(ora.up.reactome)
ora.down.BP.chrX <- match.chrX(ora.down.BP)
ora.up.BP.chrX <- match.chrX(ora.up.BP)
ora.down.hallmark.chrX <- match.chrX(ora.down.hallmark)
ora.up.hallmark.chrX <- match.chrX(ora.up.hallmark)

# check if not null or nrow > 0 and Save results
if(!is.null(ora.down.kegg.chrX) && nrow(ora.down.kegg.chrX) > 0){
    save(ora.down.kegg.chrX, file = 'ORA/ora.down.kegg.chrX.Rdata')
}
if(!is.null(ora.up.kegg.chrX) && nrow(ora.up.kegg.chrX) > 0){
    save(ora.up.kegg.chrX, file = 'ORA/ora.up.kegg.chrX.Rdata')
}
if(!is.null(ora.down.reactome.chrX) && nrow(ora.down.reactome.chrX) > 0){
    save(ora.down.reactome.chrX, file = 'ORA/ora.down.reactome.chrX.Rdata')
}
if(!is.null(ora.up.reactome.chrX) && nrow(ora.up.reactome.chrX) > 0){
    save(ora.up.reactome.chrX, file = 'ORA/ora.up.reactome.chrX.Rdata')
}
if(!is.null(ora.down.BP.chrX) && nrow(ora.down.BP.chrX) > 0){
    save(ora.down.BP.chrX, file = 'ORA/ora.down.BP.chrX.Rdata')
}
if(!is.null(ora.up.BP.chrX) && nrow(ora.up.BP.chrX) > 0){
    save(ora.up.BP.chrX, file = 'ORA/ora.up.BP.chrX.Rdata')
}
if(!is.null(ora.down.hallmark.chrX) && nrow(ora.down.hallmark.chrX) > 0){
    save(ora.down.hallmark.chrX, file = 'ORA/ora.down.hallmark.chrX.Rdata')
}
if(!is.null(ora.up.hallmark.chrX) && nrow(ora.up.hallmark.chrX) > 0){
    save(ora.up.hallmark.chrX, file = 'ORA/ora.up.hallmark.chrX.Rdata')
}