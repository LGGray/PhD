# Script to read in data from each female from OneK1k and prepare for brewerIX analysis

# Read in female IDs
female.id <- read.delim('~/datasets/OneK1k/female.ids.txt', header=F)
# Convert variantID to RefSNP rsID
vcf <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/HRC.chrX.tab')

path <- '/directflow/SCCGGroupShare/projects/lacgra/ASEReadCounter/individuals/'
result.list <- list()
for (i in female.id){
  files <- list.files(path = paste0(path, i), pattern='[ACTG]+-[0-9]+.txt', full.names = T)
  if(length(files) > 1){
    ase <- lapply(files, function(x){
      read.delim(x, colClasses = c(rep('character', 5), rep('numeric', 8)))
    })
    names(ase) <- paste0(i, ':', gsub('.txt', '', basename(files)))
    ase <- ase[sapply(ase, function(x) nrow(x) > 0)]
    ase <- lapply(ase, function(x){
      cbind(x, rsID=subset(vcf, variantID %in% x$variantID)$ID)
    })
    # Create new column: chr_pos_rs_ref_alt
    ase <- lapply(ase, function(x){
      cbind(x, chr_pos_rs_ref_alt = paste(x[,1], x[,2], x[,14], x[,4], x[,5], sep='_'))
    })
    # Transform list into dataframe
    # ase.table <- dplyr::bind_rows(ase, .id=NULL)
    # Aggregate SNP count across cells
    # ase.sum <- aggregate(x=ase.table[,-c(1,2,3,4,5,8,9,10,11,12,13,14,15)], by=list(ase.table$chr_pos_rs_ref_alt), FUN=sum)
    # Add value column
    # ase.sum <- cbind(ase.sum, value = paste(ase.sum$refCount, ase.sum$altCount, '1.0', sep=','))
    # Change column name
    # colnames(ase.sum)[1] <- 'chr_pos_rs_ref_alt'
    # Subset columns
    # ase.sum <- ase.sum[,c(1, 4)]
    # colnames(ase.sum)[2] <- i
    
    # Add count of 5 to all reads
    ase <- lapply(ase, function(x){
      x[,6:7] <- apply(x[,6:7], 2, function(x) x + 5)
      return(x)
    })
    # Create new column: value
    ase <- lapply(ase, function(x){
      cbind(x, value = paste(x[,6], x[,7], '1.0', sep=','))
    })
    # Subset for new columns
    ase.sub <- lapply(ase, function(x) x[,15:16])
    # Create one dataframe
    ase.table <- Reduce(function(df1, df2) merge(df1, df2, by = 'chr_pos_rs_ref_alt', all = TRUE), ase.sub)
    colnames(ase.table)[-1] <- names(ase.sub)
    result.list[[length(result.list)+1]] <- ase.table
    print('done')
  }
}

# Create ASER_table.txt
ase.table <- Reduce(function(df1, df2) merge(df1, df2, by = 'chr_pos_rs_ref_alt', all = TRUE), result.list)
ase.table <- apply(ase.table, 2, as.character)

# for(i in 2:ncol(ase.table)){
#   ase.table[grep('0,\\d,1.0|\\d,0,1.0',ase.table[,i]), i] <- NA
# }
# apply(ase.table[,-1], 2, function(x) x[grep('0,\\d,1.0|\\d,0,1.0',x),2])

setwd('/directflow/SCCGGroupShare/projects/lacgra/brewerIX')
write.table(ase.table, 'ASER_table.txt', row.names=F, quote = F, sep="\t")
# file.create(paste0(colnames(ase.table)[2], '_1.fq.gz'))
# file.create(paste0(colnames(ase.table)[2], '_2.fq.gz'))
# file.create(paste0(colnames(ase.table)[2], '.bam'))
# file.create(paste0(colnames(ase.table)[2], '.bam.bai'))
# file.create(paste0('do_not_run_ASER'))
#########

# setwd('/directflow/SCCGGroupShare/projects/lacgra/ASEReadCounter/individuals/684_685')
# 
# files <- list.files(pattern = '[ACTG]+-[0-9]+.txt')
# 
# ase <- lapply(files[1:100], function(x){
#   read.delim(x, colClasses = c(rep('character', 5), rep('numeric', 8)))
# })
# names(ase) <- paste0('684_685:',gsub('.txt', '', files[1:100]))
# ase <- ase[sapply(ase, function(x) nrow(x) > 0)]
# 
# 
# 
# ase <- lapply(ase, function(x){
#   cbind(x, rsID=subset(vcf, variantID %in% x$variantID)$ID)
# })
# 
# # Create new column: chr_pos_rs_ref_alt
# ase <- lapply(ase, function(x){
#   cbind(x, chr_pos_rs_ref_alt = paste(x[,1], x[,2], x[,14], x[,4], x[,5], sep='_'))
# })
# # Create new column: value
# ase <- lapply(ase, function(x){
#   cbind(x, value = paste(x[,6], x[,7], 1, sep=','))
# })
# 
# # Subset for new columns
# ase.sub <- lapply(ase, function(x) x[,15:16])
# 
# a <- dplyr::bind_rows(ase.sub)
# b <- dplyr::bind_rows(ase.sub)
# 
# ase.table <- merge(a,b, by='chr_pos_rs_ref_alt', all=T)
# colnames(ase.table) <- c('693_694', '684_685', '687_688')
# 
# # Create ASER_table.txt
# ase.table <- Reduce(function(df1, df2) merge(df1, df2, by = 'chr_pos_rs_ref_alt', all = TRUE), ase.sub)
# colnames(ase.table)[-1] <- names(ase)
# ase.table <- apply(ase.table, 2, as.character)

# dir.create('/directflow/SCCGGroupShare/projects/lacgra/test')
# setwd('/directflow/SCCGGroupShare/projects/lacgra/test')
# write.table(ase.table, 'ASER_table.txt', row.names=F, quote = F, sep="\t")
# 
# file.create(paste0(colnames(ase.table)[2], '_1.fq.gz'))
# file.create(paste0(colnames(ase.table)[2], '_2.fq.gz'))
# file.create(paste0(colnames(ase.table)[2], '.bam'))
# file.create(paste0(colnames(ase.table)[2], '.bam.bai'))
# file.create(paste0('do_not_run_ASER'))

