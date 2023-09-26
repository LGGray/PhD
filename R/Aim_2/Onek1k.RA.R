# Determine which of the OneK1k samples have RA and their respective pools

disease.meta <- readRDS('disease.meta.RDS')
metadata <- readRDS('metadata_final.Rds')

# Match 
disease.meta$X11_Rheumatoid_arthritis[is.na(disease.meta$X11_Rheumatoid_arthritis)] <- 0
RA_samples <- disease.meta$individual[disease.meta$X11_Rheumatoid_arthritis == 1]

RA_metadata <- unique(metadata[metadata$individual %in% RA_samples, c('individual', 'pool', 'sex', 'age')])
RA_metadata$pool <- as.character(RA_metadata$pool)

write.table(RA_metadata, 'RA_metadata.txt', row.names = F, quote = F)

pool_metadata <- lapply(unique(RA_metadata$pool), function(x) {
  tmp <- unique(metadata[metadata$pool == x, c('individual', 'sex', 'age')])
  tmp$RA <- ifelse(tmp$individual %in% RA_samples, TRUE, FALSE)
  return(tmp)
})
names(pool_metadata) <- unique(RA_metadata$pool)
save(pool_metadata, file = 'RA_pool_metadata.RData')



