# Function to subsample a seurat object based on another object
# df is a dataframe with cell type in first column and frequency in second

matched.sample <- function(df, seurat.object){
  seurat.list <- list()
  for(i in 1:nrow(df)){
      temp <- subset(seurat.object, predicted.celltype.l3 == df[i,1])
      temp <- subset(temp, cells = sample(Cells(temp), df[i,2]))
      seurat.list <- append(seurat.list, temp)
  }
  seurat.output <- merge(x=seurat.list[[1]], y=seurat.list[2:length(seurat.list)])
  return(seurat.output)
}



