library(pROC)
power.roc.test(auc = 0.9, null.auc = 0.5, sig.level = 0.01, 
               power = 0.90, alternative = "two.sided", ratio = 1)

power.roc.test(auc = 0.8, null.auc = 0.5, sig.level = 0.05, 
               power = 0.80, alternative = "two.sided", ratio = 1)

