# Create colour codes for cell types
library(RColorBrewer)
library(ggplot2)

# For cells
cells <- unique(gsub('.txt', '', basename(Sys.glob('*/edgeR-LRT/*.txt'))))
cellcolours <- colorRampPalette(brewer.pal(8, "Set1"))(length(cells))
# Create data frame
cell.colourdict <- data.frame(row.names=cells, colour=cellcolours)
# Create dictionary
cell.colourdict <- as.list(t(cell.colourdict)[1,])

#For studies
studies <- c('SLE', 'RA', 'UC', 'pSS', 'MS')
studycolours <- colorRampPalette(brewer.pal(5, "Set2"))(length(studies))
# Create data frame
study.colourdict <- data.frame(row.names=studies, colour=studycolours)
# Create dictionary
study.colourdict <- as.list(t(study.colourdict)[1,])

