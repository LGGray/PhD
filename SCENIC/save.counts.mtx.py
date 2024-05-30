import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('pbmc.female.control-managed.h5ad')
t=adata.raw.X.toarray()
pd.DataFrame(data=t, index=adata.obs_names, columns=adata.raw.var_names).to_csv('SCENIC/expr_matrix.csv')