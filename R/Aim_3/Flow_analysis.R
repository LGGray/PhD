library(ggplot2)
library(ggsignif)
library(ROCR)

# Read in dataframe
df <- read.delim('SLE_Lachie_Results.txt')

### Check assumptions of normality ###
model_list <- list()
for(i in (c(5:14))){
  model <- lm(df[,i] ~ df$Disease + df$Run)
  pdf(paste0('Normality_assumptions/', colnames(df)[i], '.pdf' ))
  par(mfrow=c(2,2))
  plot(model)
  dev.off()
  model_list[[colnames(df)[i]]] <- model
}

lapply(model_list, summary)

lapply(5:14, function(x){
  shapiro.test(df[,x])
})

### Plot boxplots with signifinace to show diferences in MFI and cellCount
df_wide <- reshape2::melt(df[,-c(3,4,6,12,14)])
df_wide$variable <- gsub('_Geomean', '', df_wide$variable)

pdf('markers_boxplot.pdf', width=12, height=12)
comparisons <- list(c('SLE', 'Healthy'))
ggplot(df_wide, aes(x=Disease, y=value, colour=Disease)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  theme(legend.position = 'none') +
  geom_signif(comparisons=comparisons, map_signif_level=TRUE, test='t.test', color = 'black') +
  ylab('MFI') + xlab('') +
  facet_wrap(~variable, scales = 'free_y', strip.position = "bottom")
dev.off()

pdf('cellcount_boxplot.pdf')
comparisons <- list(c('SLE', 'Healthy'))
df_cellCount <- reshape2::melt(df[,c('Disease', 'BMEM_CellCount', 'Tregs_CellCount', 'Mait_cells_Freq_Live_cells')])
ggplot(df_cellCount, aes(x=Disease, y=value, colour=Disease)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  ylab('count') + xlab('') +
  theme(legend.position = 'none') +
  geom_signif(comparisons=comparisons, map_signif_level=TRUE, test='t.test', color = 'black') +
  facet_wrap(~variable, scales = 'free_y', strip.position = "bottom") 
dev.off()

### Plot ROC curve to show predictive ability of markers ###
df$Disease <- ifelse(df$Disease == 'SLE', 1, 0)
df$Disease <- factor(df$Disease)

# PIM2
PIM2_model <- glm(Disease ~ BMEM_Geomean_PIM2 + BMEM_CellCount + Age, data = df, family = "binomial")
probabilities <- predict(PIM2_model, type = "response")
pred <- prediction(probabilities, df$Disease)
perf <- performance(pred, "tpr", "fpr")
plot(perf, col = "blue", main = "ROC for PIM2 + BMEM Cell Count + Age", xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")
auc <- performance(pred, measure = "auc")
auc_value <- auc@y.values[[1]]
# Add AUC to the plot
text(0.6, 0.2, labels = paste("AUC =", round(auc_value, 3)), col = "red", cex = 1.2)

# TIMP1
# Set up the control parameters for cross-validation
TIMP1_model <- glm(Disease ~ Tregs_Geomean_TIMP1 + Tregs_CellCount, data = df, family = "binomial")
probabilities <- predict(TIMP1_model, type = "response")
pred <- prediction(probabilities, df$Disease)
perf <- performance(pred, "tpr", "fpr")
plot(perf, col = "blue", main = "ROC for TIMP1 and Treg Cell Count", xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")
auc <- performance(pred, measure = "auc")
auc_value <- auc@y.values[[1]]
# Add AUC to the plot
text(0.6, 0.2, labels = paste("AUC =", round(auc_value, 3)), col = "red", cex = 1.2)

