library("MatrixEQTL")
library("SingleCellExperiment")
library("vcfR")
library("Seurat")
library("dplyr")

setwd("~/datasets/OneK1k/variants/MatrixEQTL/")
# read in seurat object
pbmc <- readRDS("~/datasets/OneK1k/variants/pbmc.RDS")

# Extract celltype of interest
pbmc <- subset(pbmc, predicted.celltype.l2 == "CD14 Mono")
# convert to sce object
sce <- as.SingleCellExperiment(pbmc, assay = "SCT")
# pull out sample ids
individual <- sce$individual
# pseduobulk
summed <- scuttle::aggregateAcrossCells(sce, factor(individual))

# Setting covariate file
RA <- summed$RA
age <- summed$age
pool <- summed$pool
covr <- rbind(RA, age, pool)
colnames(covr) <- colnames(summed)
write.table(covr, "covt.txt", sep = "\t", row.names = T)

# read in imputed vcf file (MAF < 0.1 & nonPAR)
vcf <- read.vcfR("~/datasets/OneK1k/variants/snpEff/all.recode.ann.vcf", verbose = F)
# extract genotype data
gt <- extract.gt(vcf)
# extract annotated data
info <- vcf@fix[,8]
# extract gene name from snpEff
genes <- sapply(strsplit(info, "\\|"), `[`, 4)

#rearrange columns to match expression data
gt <- gt[,colnames(gt) %in% colnames(summed)]
gt <- gt[,colnames(summed)]

# find and replace
for (i in 1:length(colnames(gt))){
  gt[grep("0\\|0", gt[,i]),i] <- 0
  gt[grep("1\\|0", gt[,i]),i] <- 1
  gt[grep("0\\|1", gt[,i]),i] <- 1
  gt[grep("1\\|1", gt[,i]),i] <- 2
}

#gt <- apply(gt[,1:length(colnames(gt))], 2, function(x) as.numeric(x))
gt.ann <- cbind(gt, genes)
gt.filter <- subset(gt.ann, genes %in% rownames(summed))
write.table(gt.filter[,1:length(colnames(gt.filter))-1], "snp.txt", sep = "\t", row.names = T, quote=F)

expr <- assay(summed)
expr <- subset(expr, rownames(expr) %in% gt.filter[,length(colnames(gt.filter))])
write.table(expr, "GE.txt", sep = "\t", row.names = T)

# MatrixEQTL setup
useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS
output_file_name = "RA"

pvOutputThreshold = 1e-2;
errorCovariance = numeric();

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( "snp.txt" );

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile( "GE.txt" );

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
cvrt$LoadFile( "covt.txt" );

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

# Extract lead eQTL
top_eqtls = filter(me$all$eqtls, FDR <= 0.05) %>%
  arrange(pvalue) %>%
  distinct(gene, .keep_all = TRUE)

snp_values <- read.delim("snp.txt", header=T)
colnames(snp_values) <- gsub("X", "", colnames(snp_values))
gene_values <- read.delim("GE.txt", header=T)
colnames(gene_values) <- gsub("X", "", colnames(gene_values))

# Visualise
# top_gene = "XIST"
# top_snp = "X:2835831:C:T"
top_snp = top_eqtls$snps[1]
top_gene = as.character(top_eqtls$gene[1])

top_snp_data = filter(snp_values, rownames(snp_values) == top_snp)
top_gene_data = filter(gene_values, rownames(gene_values) == top_gene)

plot_data = t(bind_rows(top_snp_data[-1], top_gene_data[-1]))
colnames(plot_data) = c("snp", "gene_expr")
plot_data = as.data.frame(plot_data)
plot_data$snp = as.factor(plot_data$snp)

lm_top = lm(plot_data[,"gene_expr"] ~ as.numeric(plot_data[,"snp"]))
summary(lm_top)

pdf(paste0(top_gene, ".eqtl.pdf"))
plot(plot_data, col="steelblue", 
     main = paste0(top_gene, " vs ", top_snp))
abline(lm_top, col="darkorange", lwd = 2, lty = 2)
y_range = range(plot_data[,"gene_expr"])
text(x=2, y=y_range[1] + 0.5*diff(y_range), paste0("p=",
                                                   format(summary(lm_top)$coefficients[2,4],
                                                          scentific=TRUE, digits=2)), col = "darkorange")
dev.off()
