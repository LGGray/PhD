library(dplyr)
library(tidyr)

pdf('cellprop.pdf')
ggplot(pbmc@meta.data, aes(x=individual, fill=cellTypist)) +
    geom_bar('position'='fill') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Total number of cells for each individual
pbmc@meta.data %>%
    group_by(individual, condition) %>%
    summarize(total_cells = n()) %>%
    ungroup() -> cell_counts

# Number of cells for each cell type
pbmc@meta.data %>%
    group_by(individual, condition, cellTypist) %>%
    summarize(celltype_counts = n()) %>%
    ungroup() -> celltype_counts
# join the cell counts and celltype counts tables
left_join(celltype_counts, cell_counts, by = c("individual", "condition")) %>%
  mutate(proportion = celltype_counts / total_cells) -> celltype_proportions
# fit a logistic regression model to compare celltype proportion between conditions
glm(proportion ~ cellTypist + condition, family = binomial(link = "logit"), data = celltype_proportions) -> logit_model
# summarize the results
summary(logit_model)