library(Seurat)
library(scProportionTest)

pbmc <- readRDS("~/datasets/OneK1k/C_vs_RA_DEout/pbmc.RDS")
data_subset <- subset(data_subset, predicted.celltype.l2 != "NA")

prop_test <- sc_utils(data_subset)

prop_test <- permutation_test(
  prop_test, cluster_identity = "predicted.celltype.l2",
  sample_1 = "1", sample_2 = "2",
  sample_identity = "sex"
)

pdf("~/datasets/OneK1k/M_vs_F_DEout/proportion.test.pdf")
permutation_plot(prop_test)
dev.off()

saveRDS(prop_test, "~/datasets/OneK1k/M_vs_F_DEout/prop_test.RDS")
