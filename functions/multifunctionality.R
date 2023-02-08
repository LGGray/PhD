# Function to calculate the multifunctionality of a gene
# Based on the method described in Gillis et al 2011 doi:10.1371/journal.pone.0017258

multifunct <- function(gene.list, gene.set){
    result <- lapply(gene.list, function(gene){
        Num_in.i <- length(grep(gene, gene.set))
        Num_out.i <- length(grep(gene, gene.set, invert = TRUE))
        return 1/(Num_in.i * Num_out.i)
    }
}