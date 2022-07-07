library(tidyverse)

setwd('/directflow/SCCGGroupShare/projects/lacgra/ASEReadCounter/individuals/693_694')

files <- list.files(pattern = '[ACTG]+-[0-9]+.txt')

ase <- lapply(files, function(x){
  read.delim(x, colClasses = c(rep('character', 5), rep('numeric', 8)))
})
names(ase) <- gsub('.txt', '', files)
ase <- ase[sapply(ase, function(x) nrow(x) > 0)]

# Create new column: chr_pos_rs_ref_alt
ase <- lapply(ase, function(x){
  cbind(x, chr_pos_rs_ref_alt = paste(x[,1], x[,2], x[,3], x[,4], x[,5], sep='_'))
})
# Create new column: value
ase <- lapply(ase, function(x){
  cbind(x, value = paste(x[,6], x[,7], 1, sep=','))
})

# Subset for new columns
ase.sub <- lapply(ase, function(x) x[,14:15])

# Create ASER_table.txt
ase.table <- Reduce(function(df1, df2) merge(df1, df2, by = 'chr_pos_rs_ref_alt', all = TRUE), ase.sub)
colnames(ase.table)[-1] <- names(ase)

write.table(ase.table, 'ASER_table.txt', row.names=F, quote = F, sep="\t")

dir.create('test')

file.create(paste0('test/', colnames(ase.table)[2], '_1.fq.gz'))
file.create(paste0('test/', colnames(ase.table)[2], '_2.fq.gz'))
file.create(paste0('test/', colnames(ase.table)[2], '.bam'))
file.create(paste0('test/', colnames(ase.table)[2], '.bai'))
file.create(paste0('test/', 'do_not_run_ASER'))



