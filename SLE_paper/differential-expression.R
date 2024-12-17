library(edgeR)
library(Seurat)


source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/Seurat2PB.R')

# Read in file from command line
pbmc <- readRDS('pbmc_female.control_managed.RDS')

for (cell in levels(pbmc)){
  # Select cell type
  print(cell)
  # subset object by cell type
  pbmc.cell <- subset(pbmc, ct_cov == cell)

  if (length(unique(pbmc.cell$condition)) < 2){
    next
  }

  # Perform Pseudobulking as per edgeR
  y <- Seurat2PB(pbmc.cell, sample='individual', cluster='condition', assay='RNA')

  # Keep samples with library size greater than 1st quartile
  # keep.samples <- y$samples$lib.size > summary(y$samples$lib.size)[2]
  # print('Number of samples dropped and retained')
  # print(table(keep.samples))
  # y <- y[, keep.samples]

  # Remove genes with low expression
  # keep.genes <- filterByExpr(y, group=y$samples$cluster)
  # print('Number of genes dropped and retained')
  # print(table(keep.genes))
  # y <- y[keep.genes, , keep=FALSE]

  # TMM normalisation
  y <- calcNormFactors(y)

  # Reorder cluster i.e condition so disease is the reference
  y$samples$cluster <- factor(y$samples$cluster, levels = c("disease", "control"))
  design <- model.matrix(~ 0 + cluster, data = y$samples)

  # Estimate gene wise dispersion
  y <- estimateDisp(y, design, robust=TRUE)

  # Fit to distribution
  fit <- glmQLFit(y, design, robust=TRUE)

  # Explicitely set up contrasts with disease as reference
  contrastMatrix <- makeContrasts(DiseaseVsControl = clusterdisease - clustercontrol, levels = design)
  
  # Perform QLFTest
  qlf <- glmQLFTest(fit, contrast = contrastMatrix)
  print(summary(decideTests(qlf)))

  # Extract results
  res = topTags(qlf, n = Inf)[[1]]

  # Save to file
  cell = gsub("/|-| ", "_", cell)
  write.table(res, paste0("differential.expression/edgeR/", cell, ".txt"),
              row.names=F, sep="\t", quote = F)
}

print("Done with edgeR-QLF")

