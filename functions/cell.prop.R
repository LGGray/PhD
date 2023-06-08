library(dplyr)
library(tidyr)
library(tidyverse)
library(emmeans)
library(ggplot2)
library(ggrepel)
library(Seurat)

pbmc <- readRDS("pbmc.female.RDS")

metadata <- pbmc@meta.data

cell_type_counts <- metadata[,c("condition", "individual", "cellTypist")]

cell_type_counts %>% count(individual) -> total_counts

left_join(cell_type_counts, total_counts, 
          "individual") -> cell_type_counts

cell_type_counts %>% 
  group_by(condition, individual, n) %>% 
  count(cellTypist, name='count') -> cell_type_counts

cell_type_counts$other <- cell_type_counts$n - cell_type_counts$count

formula = cbind(count, other) ~ cellTypist * condition
model <- glm(formula = formula, family = 'binomial', data = cell_type_counts)

summary(model)

anova(model, test = "Chisq")

emm <- emmeans(model, specs = ~ condition | cellTypist)
con <- contrast(emm, method="revpairwise")
summary(con)

df <- data.frame(con)
df$estimate <- round(df$estimate, 2)
df[,-c(1,4,5)]

write.table(data.frame(con), 'cell.proportions.txt', sep='\t', quote=FALSE, row.names=FALSE)

plot.data <- data.frame(con)
plot.data <- plot.data[!is.na(plot.data$p.value),]
pdf('cell.proportions.pdf')
ggplot(plot.data, aes(x=estimate, y=-log10(p.value))) +
  geom_point() +
  geom_text_repel(aes(label=cellTypist), size=5, color = 'black') +
  theme_minimal() +
  labs(y = '-log10(p-value)', title = 'Cell type proportion change', x = 'estimate')
dev.off()

pdf('cellprop.pdf')
ggplot(pbmc@meta.data, aes(x=individual, fill=cellTypist)) +
    geom_bar('position'='fill') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#### edgeR method
library(edgeR)
library(Seurat)
library(MAST)
library(qvalue)
library(tidyverse)
library(ggplot2)
library(ggrepel)

pbmc <- readRDS('pbmc.female.RDS')

abundances <- table(pbmc$cellTypist, pbmc$individual) 
abundances <- unclass(abundances) 
head(abundances)

# Attaching some column metadata.
extra.info <- pbmc@meta.data[match(colnames(abundances), pbmc$individual),]
y.ab <- DGEList(abundances, samples=extra.info[,"condition"])
y.ab

# Filter out low-abundance labels
keep <- filterByExpr(y.ab, group=y.ab$samples$samples)
y.ab <- y.ab[keep,]
summary(keep)

design <- model.matrix(~0+samples, y.ab$samples)
contrasts <- makeContrasts(disease_vs_control = samplesdisease - samplescontrol, levels = design)
contrast_matrix <- contrasts[ ,c("disease_vs_control")]
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)

res <- glmQLFTest(fit.ab, contrast=contrast_matrix)
summary(decideTests(res))
topTags(res)

res = topTags(res, n = Inf) %>%
    as.data.frame() %>%
    # mutate(FDR=qvalue(p = PValue)$qvalues) %>%
    rownames_to_column('celltype')
    
# Plot volcano plot
res$threshold <- ifelse(res$FDR < 0.05 & abs(res$logFC) > 1, TRUE, FALSE)
pdf('differential.abundance.volcano.pdf')
ggplot(res, aes(x=logFC, y=-log10(FDR), colour=threshold)) +
    geom_point() +
    geom_label_repel(aes(label=celltype), size=3, nudge_x=0.5) +
    theme_bw() +
    labs(x = "log2 fold change", y = "-log10 FDR")
dev.off()