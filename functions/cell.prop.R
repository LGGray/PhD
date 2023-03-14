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