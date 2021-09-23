library(ggplot2)
library(reshape2)
library("dplyr")
library(ggdendro)
library(grid)
library(MASS)

setwd("~/external/ClusterHome/datasets/OneK1k/variants/MatrixEQTL/")

# Read in files and fox colnames
snp_values <- read.delim("snp.txt", header=T)
colnames(snp_values) <- gsub("X", "", colnames(snp_values))
gene_values <- read.delim("GE.txt", header=T)
colnames(gene_values) <- gsub("X", "", colnames(gene_values))
covt <- read.delim("covt.txt", header=T)
colnames(covt) <- gsub("X", "", colnames(covt))
eqtls <- read.delim("RA.txt")

# Filter for significant eQTLs
top_eqtls = filter(eqtls, FDR <= 0.05) %>%
  arrange(FDR) %>%
  distinct(gene, .keep_all = TRUE)

# Order on RA status 0 or 1
covt.order <- covt[,order(t(covt)[,3])]
snp.sig <- subset(snp_values, rownames(snp_values) %in% top_eqtls$SNP)
snp.sig <- snp.sig[,colnames(covt.order)]
snp.sig$id <- rownames(snp.sig)

# Generate dendrogram for heatmap
snp.dendo <- as.dendrogram(hclust(d = dist(x = snp.sig)))
dendro.plot <- ggdendrogram(data = snp.dendo, rotate = TRUE)
dendro.plot <- dendro.plot + theme(axis.text.y = element_text(size = 0))
dendro.plot

# Reorder genotype data so it matches dendrogram
snp.sig.order <- order.dendrogram(snp.dendo)
snp.sig.long$snp <- factor(x = snp.sig.long$snp,
                           levels = rownames(snp.sig)[snp.sig.order], 
                           ordered = TRUE)
# Long format melt
snp.sig.long <- melt(snp.sig, value.name="id")
colnames(snp.sig.long) <- c("snp", "sample", "genotype")

# Generate heatmap 
heatmap.plot <- ggplot(snp.sig.long, aes(x=sample, y=snp, fill=genotype)) + geom_tile(colour="white") +
  scale_fill_gradient2(low = "purple", high = "red", mid = "blue", 
                       midpoint = 1, limit = c(0,2), space = "Lab", 
                       name="Genotype") +
  ggtitle("eQTL genotype enrichment") +
  xlab("Samples") + ylab("") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "top")

# Add heatmap together and dendrogram together. 
grid.newpage()
print(heatmap.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.80, y = 0.44, width = 0.2, height = 0.93))

# Filter gene expression for eQTL genes
gene.sig <- subset(gene_values, rownames(gene_values) %in% top_eqtls$gene)
gene.sig$gene <- rownames(gene.sig)
# Melt to long format
gene.sig.long <- melt(gene.sig, id = "gene")

# Plot boxplot 
ggplot(gene.sig.long, aes(x=gene, y=value, fill=gene)) + geom_boxplot(varwidth = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 6, hjust = 1)) +
  theme(legend.position = "none") +
  ggtitle("Top eQTL gene expression") +
  xlab("Genes") + ylab("Psudobulked expression")


# Plot gene expression vs genotype
target = subset(top_eqtls, gene == "PRDX4")
top_snp = target$SNP
top_gene = as.character(target$gene)
RA_cond = covt[3,]

top_snp_data = filter(snp_values, rownames(snp_values) == top_snp)
top_gene_data = filter(gene_values, rownames(gene_values) == top_gene)

plot_data = t(bind_rows(top_snp_data, top_gene_data, RA_cond))
colnames(plot_data) = c("snp", "gene_expr", "RA")
plot_data = as.data.frame(plot_data)
plot_data$snp = as.factor(plot_data$snp)
plot_data$RA = as.factor(plot_data$RA)

# Perform linear regression
lm_top = lm(plot_data[,"gene_expr"] ~ as.numeric(plot_data[,"snp"]) + plot_data[,"RA"])
summary(lm_top)

# GLM with negative binomial 
#summary(m1 <- glm.nb(plot_data[,"gene_expr"] ~ as.numeric(plot_data[,"snp"]) + plot_data[,"RA"]))

#pdf(paste0(top_gene, ".eqtl.pdf"))
plot(plot_data[,1:2], col="steelblue", 
     main = paste0(top_gene, " vs ", top_snp))
abline(lm_top, col="darkorange", lwd = 2, lty = 2)
y_range = range(plot_data[,"gene_expr"])
text(x=2, y=y_range[1] + 0.5*diff(y_range), paste0("p=",
                                                   format(summary(lm_top)$coefficients[2,4],
                                                          scentific=TRUE, digits=2)), col = "darkorange")
#dev.off()

