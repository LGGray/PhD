# Take command line arguments
study <- commandArgs(trailingOnly=TRUE)

# Read in files
metrics1 <- read.delim(paste0(study[1], '/exp.matrix/metrics/Metrics.combined.txt'))
metrics2 <- read.delim(paste0(study[2], '/exp.matrix/metrics/Metrics.combined.txt'))
# combine files
merged <- merge(metrics1, metrics2, by='model')

# Subset for chrX models
merged <- merged[grep('.chrX', merged$model),]

# Order by largest F1 score for both models
merged <- merged[order(merged$F1.x, merged$F1.y, decreasing=TRUE),]
# Subset for F1 score > 0.8
merged <- merged[merged$F1.x > 0.8 & merged$F1.y > 0.8,]

write.table(merged, file=paste('analysis', study[1], study[2], 'best.model.txt', sep='.'), sep='\t', row.names=FALSE, quote=FALSE)