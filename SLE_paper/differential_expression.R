library(edgeR)
library(Seurat)

source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/Seurat2PB.R')

# Read in Seurat file
pbmc <- readRDS('pbmc_female.control_managed.RDS')

for (cell in levels(pbmc)){
  # Select cell type
  print(cell)
  # subset object by cell type
  pbmc.cell <- subset(pbmc, cell_type_detailed == cell)

  if (length(unique(pbmc.cell$condition)) < 2){
    next
  }

  # Perform Pseudobulking as per edgeR
  y <- Seurat2PB(pbmc.cell, sample='ind_cov', cluster='disease', assay='COMBAT_LogNorm')

  # Filter out lowly expressed genes
  keep.genes <- filterByExpr(y, group=y$samples$cluster)
  y <- y[keep.genes, , keep=FALSE]

  # TMM normalisation
  y <- calcNormFactors(y)

  # Reorder cluster i.e condition so disease is the reference
  y$samples$cluster <- factor(y$samples$cluster, levels = c("disease", "control"))

  batch <- factor(y$samples$batch)
  cluster <- factor(y$samples$cluster)
  design <- model.matrix(~ + batch + cluster)

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
  cell = gsub("+|-| ", "_", cell)
  write.table(res, paste0("edgeR/", cell, ".txt"),
              row.names=F, sep="\t", quote = F)
}

print("Done with edgeR-QLF")

