library(karyoploteR)
library(GenomicRanges)

# Load a GTF file (downloaded from Ensembl)
gtf <- rtracklayer::import("~/Downloads/Homo_sapiens.GRCh38.112.gtf.gz")

# Load features and background genes
features <- read.csv('~/Dropbox (Garvan)/LGray PhD/Thesis Results/Aim_2/top_chrX.consistent.csv')
load('~/Dropbox (Garvan)/LGray PhD/Thesis Results/Aim_2/selected_features.Rdata')
features <- unique(unlist(selected_features$all_features[grepl('.chrX', names(selected_features$all_features))]))

background <- read.delim('~/Dropbox (Garvan)/LGray PhD/Thesis Results/chrX_biomaRt.txt')
background <- subset(background, Gene.name != '')
background <- background$Gene.name[!background$Gene.name %in% features]

# Extract coordinates for features and background
feature_coords <- subset(gtf, gene_name %in% features & type == 'gene')
feature_coords_df <- as.data.frame(feature_coords)
feature_coords_df$seqnames <- paste0('chr', feature_coords_df$seqnames)
feature_coords_gr <- GRanges(seqnames = feature_coords_df$seqnames,
                             ranges = IRanges(
                               start = feature_coords_df$start,
                               end = feature_coords_df$end
                             ))

background_coords <- subset(gtf, gene_name %in% background & type == 'gene')
background_coords_df <- as.data.frame(background_coords)
background_coords_df$seqnames <- paste0('chr', background_coords_df$seqnames)
background_coords_gr <- GRanges(seqnames = background_coords_df$seqnames,
                                ranges = IRanges(
                                  start = background_coords_df$start,
                                  end = background_coords_df$end
                                ))

# Define XIST locus
xist_locus <- subset(gtf, gene_name == "XIST" & type == 'gene')
xist_start <- start(xist_locus)

line <- read.delim('~/Downloads/rmsk.txt', header=F)
line <- subset(line, V6 == 'chrX' & V12 == 'LINE')
line <- GRanges(seqnames = line$V6,
                ranges = IRanges(
                  start = line$V7,
                  end = line$V8
                ))

# Calculate distances from XIST
calculate_distance <- function(gene_coords, target_coords) {
  sapply(seq_along(gene_coords), function(i) {
    min(abs(start(gene_coords[i]) - start(target_coords)), 
        abs(end(gene_coords[i]) - end(target_coords)))
  })
}

# Calculate the distance to the closest LINE element
calculate_distance_to_closest_line <- function(gene_coords, line_coords) {
  sapply(seq_along(gene_coords), function(i) {
    min_distance <- min(abs(start(gene_coords[i]) - start(line_coords)),
                        abs(start(gene_coords[i]) - end(line_coords)),
                        abs(end(gene_coords[i]) - start(line_coords)),
                        abs(end(gene_coords[i]) - end(line_coords)))
    return(min_distance)
  })
}

feature_to_xist <- calculate_distance(feature_coords_gr, xist_locus)
feature_to_closest_line <- calculate_distance_to_closest_line(feature_coords_gr, line)

background_to_xist <- calculate_distance(background_coords_gr, xist_locus)
background_to_closest_line <- calculate_distance_to_closest_line(background_coords_gr, line)

set.seed(123)
# Function to perform a permutation test
permutation_test <- function(feature_distances, background_distances, num_permutations = 10000) {
  observed_diff <- median(feature_distances) - median(background_distances)
  
  all_distances <- c(feature_distances, background_distances)
  n_features <- length(feature_distances)
  
  perm_diffs <- replicate(num_permutations, {
    perm_sample <- sample(all_distances, n_features)
    perm_diff <- median(perm_sample) - median(setdiff(all_distances, perm_sample))
    return(perm_diff)
  })
  
  p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
  
  return(p_value)
}

# Calculate p-values using permutation test
xist_permutation_p_value <- permutation_test(feature_to_xist, background_to_xist)
line_permutation_p_value <- permutation_test(feature_to_line, background_to_line)








cpg <- read.delim('~/Downloads/cpgIslandExt.txt', header=F)
cpg <- subset(cpg, V2 == 'chrX')
cpg <- GRanges(seqnames = cpg$V2,
               ranges = IRanges(
                 start = cpg$V3,
                 end = cpg$V4
               ))


gene_density <- as.data.frame(subset(gtf, seqnames=='X', type=='gene'))
gene_density$seqnames <- paste0('chr', gene_density$seqnames)
gene_density <- GRanges(seqnames = gene_density$seqnames,
                        ranges = IRanges(
                          start = gene_density$start,
                          end = gene_density$end
                        ))

kp <- plotKaryotype(plot.type=2, chromosomes = "chrX")
kpPlotDensity(kp, data=line, r0=0, r1=0.3)
#kpPlotDensity(kp, data=cpg, data.panel = 2, r0=0, r1=0.3)
kpPlotMarkers(kp, chr=feature_coords_df$seqnames, x=feature_coords_df$start, labels=feature_coords_df$gene_name)
kpPlotDensity(kp, data=gene_density, data.panel=2)

kp <- plotKaryotype(genome=subset(gtf, seqnames == 'X'))

kpPlotRegions(kp, line[,c('V6', 'V7', 'V8')], col="#AACCFF", r0=0, r1=0.25)
kpPlotRegions(kp, cpg[,c('V2', 'V3', 'V4')], col="#FFCCAA", r0=0.3, r1=0.55)
