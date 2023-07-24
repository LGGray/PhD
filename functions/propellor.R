library(speckle)
library(limma)
library(ggplot2)
library(Seurat)

pbmc <- readRDS(commandArgs(trailingOnly = TRUE)[1])
disease <- commandArgs(trailingOnly = TRUE)[2]

p1 <- plotCellTypeProps(clusters = pbmc$cellTypist, sample = pbmc$individual) + theme(axis.text.x = element_text(angle = 45))+ ggtitle("Refined cell type proportions") + 
theme(plot.title = element_text(size = 18, hjust = 0))
p1 + theme_bw() + theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust=1))                       
pdf('pbmc_celltype_props.pdf', width=10, height=10)
p1
dev.off()

#output.logit <- propeller(clusters=pbmc$cellTypist, sample=pbmc$individual, group=pbmc$condition, transform='logit')
output.asin <- propeller(clusters=pbmc$cellTypist, sample=pbmc$individual, group=pbmc$condition, transform='asin')

# Plot scatterplot of T-statistic vs -log10(FDR). colour and label if significant
pdf('pbmc_celltype_props_diff_abundance.pdf')
ggplot(output.asin, aes(x=Tstatistic, y=-log10(FDR))) + 
    geom_point(aes(colour=ifelse(FDR<0.05, "red", "black"))) +
    geom_text(aes(label=ifelse(FDR<0.05, rownames(output.asin), '')), position = position_jitter(width = 1.5, height = 1.5, ), vjust="inward",hjust="inward") +
    ylab("-log10(FDR)") + xlab("T-statistic") + 
    # rename the legend
    scale_colour_manual(name="Significant", values=c("black", "red"), labels=c("No", "Yes")) +
    ggtitle(paste(disease, " Differential Abundance"))
dev.off()

saveRDS(output.asin, 'propellor.asin.RDS') 


