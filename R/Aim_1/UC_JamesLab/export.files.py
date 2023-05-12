# import scanpy as sc
# import pandas as pd

# # Load the data
# adata = sc.read_h5ad('uc_healthy_only_jameslab_og.h5ad')
# # adata = adata.raw.to_adata()

# metadata = pd.DataFrame(adata.obs)
# metadata.to_csv('metadata.csv')
# t=adata.X.toarray()
# pd.DataFrame(data=t, index=adata.obs_names, columns=adata.var_names).to_csv('exp_counts.csv')

import scanpy as sc
from scipy.sparse import save_npz
import pandas as pd

# Load the data
adata = sc.read_h5ad('uc_healthy_only_jameslab_og.h5ad')
# Export sparse matrix
save_npz('exp_counts.npz', adata.X)
metadata = pd.DataFrame(adata.obs)
metadata.to_csv('metadata.tsv', sep='\t', index=True)
metadata['batch'] = metadata['svs']
metadata.to_csv('cell_batch.tsv', sep='\t', index=True)
features = pd.DataFrame(adata.var)
features['gene_symbol'] = features.index
features = features[['gene_ids', 'gene_symbol', 'feature_types']]
features.to_csv('features.tsv', sep='\t', index=False)
barcodes = pd.DataFrame(adata.obs.index)
barcodes.to_csv('barcodes.tsv', sep='\t', index=False, header=False)