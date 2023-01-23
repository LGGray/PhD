library(tidyverse)
library(ggrepel)
library(emmeans)
library(ggplot2)
library(ggrepel)
library(Seurat)

setwd('~/datasets/lupus.Chun')

args = commandArgs(trailingOnly=TRUE)
args <- as.numeric(args)

lst <- c('european', 'asian')

pbmc <- readRDS(paste0('seurat.object/pbmc.', lst[args], '.RDS'))

metadata <- pbmc@meta.data

metadata$development_stage <- gsub('-year-old human stage', '', 
                                   metadata$development_stage)
metadata$development_stage <- as.numeric(metadata$development_stage)
colnames(metadata)[34] <- 'age'
metadata$disease <- gsub('systemic lupus erythematosus', 'SLE', 
                         metadata$disease)

cell_type_counts <- metadata[,c("disease", "age", "sample_uuid", 
                                "predicted.celltype.l2")]

cell_type_counts %>% count(sample_uuid) -> total_counts

left_join(cell_type_counts, total_counts, 
          "sample_uuid") -> cell_type_counts
colnames(cell_type_counts)[3] <- "individual"
colnames(cell_type_counts)[4] <- "cluster"
colnames(cell_type_counts)[5] <- "total"

cell_type_counts %>% group_by(disease, age, individual, total) %>% 
  count(cluster) -> cell_type_counts

cell_type_counts <- cell_type_counts[,c(1, 2, 3, 5, 4, 6)]

cell_type_counts$other <- cell_type_counts$total - cell_type_counts$n

cell_type_counts %>% filter(total > 0) -> df_a

model0 <- glm(
  formula = cbind(n, other) ~ cluster,
  family = binomial(link = 'logit'),
  data = df_a
)

emm0 <- emmeans(model0, specs = ~ cluster)
emm0 %>%
  summary(infer = TRUE, type = 'response') %>%
  arrange(prob) -> cell_type_probs

df_a %>% filter(disease %in% c('normal', 'SLE')) -> df
formula = cbind(n, other) ~ cluster * disease + cluster 
model1 <- glm(formula = formula, family = 'binomial', data = df)

summary(model1)

emm1 <- emmeans(model1, specs = revpairwise ~ disease | cluster)
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

write.table(m_results, paste0('cell.prop.GLM.', lst[args],'.txt'), row.names=F, quote=F, sep='\t')

options(ggrepel.max.overlaps = Inf)
pdf(paste0('cell.prop.GLM.', lst[args],'.pdf'))
(
  ggplot(aes(x = prob, y = odds.ratio, color = p.value < 0.05), data = m_results)
  + geom_point()
  + geom_text_repel(aes(label = cluster), color = 'black', data = m_results %>% filter(abs(log(odds.ratio)) > log(1)))
  + scale_x_log10()
  + scale_y_log10()
  + theme_minimal()
  + labs(y = 'SLE / normal (odds ratio)', title = 'Cell type proportion change', x = 'Average abundance (probability)')
)
dev.off()
