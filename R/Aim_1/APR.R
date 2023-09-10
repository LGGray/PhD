library(Seurat)
library(speckle)
library(ggplot2)
library(reshape2)

pbmc <- readRDS('pbmc.female.RDS')

# Perform propellor cell type abundance testing
output.asin <- propeller(clusters=pbmc$cellTypist, sample=pbmc$individual, group=pbmc$condition, transform='asin')
colnames(output.asin)[1] <- 'celltype'

# Calculate percentage of each celltype in individual 
cellCount <- as.data.frame.matrix(table(pbmc$cellTypist, pbmc$individual))
cellperc <- cellCount / colSums(cellCount) * 100
cellperc <- cbind(celltype=rownames(cellperc), cellperc)

cellperc.melt <- melt(cellperc)
cellperc.melt$condition <- ifelse(grepl('HC', cellperc.melt$variable), 'control', 'disease')

# Plot percentage of each celltype in each individual
cell_order <- c(
  "CRTAM+ gamma-delta T cells", "MAIT cells", "Regulatory T cells", 
  "Tcm/Naive cytotoxic T cells", "Tcm/Naive helper T cells", 
  "Tem/Effector helper T cells", "Tem/Temra cytotoxic T cells", 
  "Tem/Trm cytotoxic T cells", "Naive B cells", "Memory B cells", 
  "Plasma cells", "CD16+ NK cells", "NK cells", "Classical monocytes", 
  "Non-classical monocytes", "DC2", "pDC", "Late erythroid", 
  "Megakaryocyte precursor", "Megakaryocytes/platelets", "Mast cells"
)
cellperc.melt$celltype <- factor(cellperc.melt$celltype, levels=cell_order)
pdf('APR/cellperc.boxplot.pdf', width=10, height=10)
ggplot(cellperc.melt, aes(x=celltype, y=value, fill=condition)) + 
    geom_boxplot() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = unit(c(1, 1, 1, 3), "lines")) + # Adjust the right margin
    labs(x='', y='Percentage of cells (%)', fill='Condition')
dev.off()

# Type I Interferon
ISG <- c(
    "ISG15", "IFI6", "IFI44L", "IFI44", "RSAD2", "CXCL10", "IFIT2", "IFIT3", "IFIT1", "IFITM3", "OAS1", "OAS3", "OAS2", "OASL", "EPSTI1", 
    "RNASE1", "RNASE2", "IFI27", "XAF1", "LGALS3BP", "SIGLEC1", "USP18", "APOBEC3A", "APOBEC3B", "MX1"
)

