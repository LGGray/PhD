library(Seurat)
library(progeny)
library(dplyr)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)

pbmc <- readRDS("/directflow/SCCGGroupShare/projects/lacgra/seurat.object/pbmc.all.RDS")
pbmc.subset <- subset(pbmc, features = VariableFeatures(pbmc), subset = RA == "Y")
# Compute progeny activity score and add to suerat object as progeny assay
pbmc.subset <- progeny(pbmc.subset, scale=F, organism="Human", top=500, perm=1, assay_name = "SCT",
                       return_assay=T)
# Scale pathway activity scores
pbmc.subset <- Seurat::ScaleData(pbmc.subset, assay = "progeny")
# Transform progeny scores into a dataframe
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(pbmc.subset, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)
# Cells which belong to each cluster
CellsClusters <- data.frame(Cell = names(Idents(pbmc.subset)), 
                            CellType = as.character(Idents(pbmc.subset)),
                            stringsAsFactors = FALSE)
# Match progeny score with the cell clusters
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)
# summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

pdf("~/plots/OneK1K/progeny/RA.progeny.pdf")

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))

progeny_hmap = pheatmap(t(summarized_progeny_scores_df),fontsize=8, 
                        fontsize_row = 8, 
                        color=myColor, breaks = progenyBreaks,
                        cluster_rows=F, cluster_cols=F,
                        main = "Rheumatoid arthritis", angle_col = 90,
                        treeheight_col = 0,  border_color = NA)
dev.off()

HC.data <- summarized_progeny_scores_df
RA.data <- summarized_progeny_scores_df

diff.data <- HC.data - RA.data
diff.data <- log(diff.data+1)


pdf("~/plots/OneK1K/progeny/difference.progeny.pdf")

progenyBreaks = c(seq(min(diff.data), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(diff.data)/paletteLength, 
                      max(diff.data), 
                      length.out=floor(paletteLength/2)))

progeny_hmap = pheatmap(t(diff.data),fontsize=8, 
                        fontsize_row = 8, 
                        color=myColor, breaks = progenyBreaks,
                        cluster_rows=F, cluster_cols=F,
                        main = "Control - Rheumatoid arthritis", angle_col = 90,
                        treeheight_col = 0,  border_color = NA)
dev.off()
