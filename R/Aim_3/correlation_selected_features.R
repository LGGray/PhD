library(dplyr)

## Testing for correlation between selected features
split_list <- list()
for (split in 1:10){

    feature.files <- list.files(paste0('split_', split, '/features'), pattern='intersection', full.names=TRUE)
    features <- lapply(feature.files, function(x) read.csv(x)$Feature)
    names(features) <- gsub('intersection_|.csv', '', basename(feature.files))

    features <- features[sapply(features, function(x) length(x) > 1)]

    correlation_list <- list()
    for(i in 1:length(features)) {
        mtx <- readRDS(paste0(names(features)[i], '.RDS'))
        mtx <- mtx[, features[[i]]]
        # drop ancestry from mtx
        mtx <- mtx[,colnames(mtx) != 'ancestry']
        correlation <- cor(mtx, method='spearman')
        tmp <- data.frame(
            celltype=names(features)[i],
            min=min(correlation[upper.tri(correlation)]), 
            max=max(correlation[upper.tri(correlation)]), 
            mean=mean(correlation[upper.tri(correlation)])
        )
        correlation_list[[names(features)[i]]] <- tmp
    }
    names(correlation_list) <- names(features)
    correlation_list <- dplyr::bind_rows(correlation_list, .id=NULL)
    split_list[[split]] <- correlation_list
}

correlation_df <- dplyr::bind_rows(split_list, .id='split')
write.csv(correlation_list, 'figures/correlation_df.csv')

correlation_df %>% 
    group_by(celltype) %>% 
    summarise(
        min=mean(min), 
        max=mean(max), 
        mean=mean(mean)
    ) %>%
    arrange(desc(max))
