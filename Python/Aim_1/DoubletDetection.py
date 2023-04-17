import scanpy as sc
import doubletdetection
import numpy as np
import pandas as pd

adata = sc.read_10x_mtx('.')
barcodes = pd.read_csv('barcodes.tsv.gz', header=None)

# remove "empty" genes
sc.pp.filter_genes(adata, min_cells=1)

clf = doubletdetection.BoostClassifier(
    n_iters=50,
    clustering_algorithm="louvain",
    standard_scaling=True,
    pseudocount=0.1,
    n_jobs=-1,
)
doublets = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
doublet_score = clf.doublet_score()

adata.obs["doublet"] = doublets
adata.obs["doublet_score"] = doublet_score

# Convergence of coublet calls
f = doubletdetection.plot.convergence(clf, save='convergence_test.pdf', show=True, p_thresh=1e-16, voter_thresh=0.5)
#Doublets on UMAP
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
# save image
sc.pl.umap(adata, color=["doublet", "doublet_score"], save='doublet_detection.pdf')
sc.pl.violin(adata, "doublet_score", save='doublet_score.pdf')

results = pd.Series(doublets, name="DoubletDetection_DropletType")
dataframe = pd.concat([barcodes, results], axis=1)
dataframe.columns = ["CellID", "DoubletDetection_DropletType"]
dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(1.0, "doublet")
dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(0.0, "singlet")
dataframe.to_csv('DoubletDetection_doublets_singlets.tsv', sep = "\t", index = False)
