library(tidyverse)
library(ggrepel)
library(emmeans)
library(ggplot2)
library(ggrepel)

pbmc <- readRDS('datasets/GSE193770/pbmc.female.RDS')

metadata <- pbmc@meta.data

cell_type_counts <- metadata[,c("condition", "sample", "predicted.celltype.l2")]

cell_type_counts %>% count(sample) -> total_counts

left_join(cell_type_counts, total_counts, 
          "sample") -> cell_type_counts
colnames(cell_type_counts)[3] <- "cluster"
colnames(cell_type_counts)[4] <- "total"

cell_type_counts %>% group_by(condition, sample, total) %>% 
  count(cluster) -> cell_type_counts
colnames(cell_type_counts)[5] <- 'count'

cell_type_counts$other <- cell_type_counts$total - cell_type_counts$count

cell_type_counts <- cell_type_counts[,c(1, 2, 4, 5, 3, 6)]

cell_type_counts %>% filter(total > 0) -> df_a

model0 <- glm(
  formula = cbind(count, other) ~ cluster,
  family = binomial(link = 'logit'),
  data = df_a
)

emm0 <- emmeans(model0, specs = ~ cluster)
emm0 %>%
  summary(infer = TRUE, type = 'response') %>%
  arrange(prob) -> cell_type_probs

df_a %>% filter(condition %in% c('MS', 'HC')) -> df
formula = cbind(count, other) ~ cluster * condition
model1 <- glm(formula = formula, family = 'binomial', data = df)

summary(model1)

emm1 <- emmeans(model1, specs = revpairwise ~ condition | cluster)
emm1$contrasts %>%
  summary(infer = TRUE, type = 'response') %>%
  rbind() %>%
  as.data.frame() -> c_results

c_results %>% arrange(desc(odds.ratio)) -> celltypeGLM.OR

emm2 <- emmeans(model1, specs = ~ cluster)
emm2 %>%
  summary(type = 'response') %>%
  select(cluster, prob) -> mean_probs

c_results %>% left_join(mean_probs) -> m_results

pdf("datasets/GSE193770/cellprop.GLM.pdf")
(
  ggplot(aes(x = prob, y = odds.ratio, color = p.value < 0.05), data = m_results)
  + geom_point()
  + geom_text_repel(aes(label = cluster), size=5, color = 'black', data = m_results %>% filter(abs(log(odds.ratio)) > log(1)))
  + scale_x_log10()
  + scale_y_log10()
  + theme_minimal()
  + labs(y = 'MS / HC (odds ratio)', title = 'Cell type proportion change', x = 'Average abundance (probability)')
)
dev.off()
