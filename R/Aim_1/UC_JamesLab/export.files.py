import scanpy as sc
import pandas as pd

# Load the data
adata = sc.read_h5ad('uc_vs_healthy_jameslab.h5ad')
adata = adata.raw.to_adata()

metadata = pd.DataFrame(adata.obs)
metadata.to_csv('metadata.csv')
counts = pd.DataFrame(adata.X.todense(), columns=adata.var.index.values, index=adata.obs.index.values)
counts.transpose().to_csv('exp_counts.csv')