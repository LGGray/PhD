library(data.table)
library(magrittr)
library(Seurat)


# save metadata as datatable
md <- data_subset@meta.data %>% as.data.table

cellcounts <- md[, .N, by = c("sex", "predicted.celltype.l2")]

cellcounts2 <- md[, .N, by = c("sex", "predicted.celltype.l2")] %>% 
  dcast(., sex ~ predicted.celltype.l2, value.var = "N")

write.table(data, "external/ClusterHome/datasets/OneK1k/M_vs_F_DEout/cellcounts2.txt", 
            sep="\t", row.names=F)
