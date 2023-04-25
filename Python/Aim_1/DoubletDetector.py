import scanpy as sc
import doubletdetection
import numpy as np
import pandas as pd

# Read in expression matrix
adata = sc.read_10x_mtx('.')
# Add metadata file
metadata = pd.read_csv('cell_batch.tsv.gz', sep='\t', header=0)
adata.obs['individual'] = metadata['batch'].values
# Read in barcodes file
barcodes = pd.read_csv('barcodes.tsv.gz', header=None)

# remove "empty" genes
sc.pp.filter_genes(adata, min_cells=1)

# run doubletdetection over each individual sample
adata_list = []
for individual in adata.obs['individual'].unique():
    adata_sub = adata[adata.obs['individual'] == individual]
    print('processing individual', individual)
    clf = doubletdetection.BoostClassifier(n_iters=50, clustering_algorithm="louvain", standard_scaling=True, pseudocount=0.1, n_jobs=-1,)
    
    doublets = clf.fit(adata_sub.X).predict(p_thresh=1e-16, voter_thresh=0.5)
    doublet_score = clf.doublet_score()

    f = doubletdetection.plot.convergence(clf, save='convergence_test'+individual+'.pdf', show=True, p_thresh=1e-16, voter_thresh=0.5)


    adata_sub.obs["doublet"] = doublets
    adata_sub.obs["doublet_score"] = doublet_score

    adata_list.append(adata_sub)

adata = adata_list[0].concatenate(adata_list[1:], )

barcodes = pd.DataFrame(adata.obs.index, columns=['CellID'])

results = pd.Series(adata.obs["doublet"], name="DoubletDetection_DropletType")
barcodes = barcodes.reset_index(drop=True)
results = results.reset_index(drop=True)
dataframe = pd.concat([barcodes['CellID'], results], axis=1)
dataframe.columns = ["CellID", "DoubletDetection_DropletType"]
dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(1.0, "doublet")
dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(0.0, "singlet")
dataframe.to_csv('DoubletDetection_doublets_singlets.tsv', sep = "\t", index = False)