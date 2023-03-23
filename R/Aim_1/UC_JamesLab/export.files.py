import scanpy as sc
import pandas as pd

# Load the data
adata = sc.read_h5ad('uc_healthy_only_jameslab_og.h5ad')
# adata = adata.raw.to_adata()

metadata = pd.DataFrame(adata.obs)
metadata.to_csv('metadata.csv')
t=adata.X.toarray()
pd.DataFrame(data=t, index=adata.obs_names, columns=adata.var_names).to_csv('exp_counts.csv')