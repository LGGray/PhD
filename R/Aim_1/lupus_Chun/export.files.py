import scanpy as sc
from scipy.sparse import save_npz
import pandas as pd

# Load the data
adata = sc.read_h5ad('local.h5ad')
# Export sparse matrix
save_npz('exp_counts.npz', adata.raw.X)
metadata = pd.DataFrame(adata.obs)
metadata['batch'] = metadata['sample_uuid']
metadata.to_csv('cell_batch.tsv', sep='\t', index=True)
features = pd.DataFrame(adata.var)
features['gene_id'] = features.index
features = features[['gene_id', 'feature_name', 'feature_biotype']]
features.to_csv('features.tsv', sep='\t', index=False, header=False)
barcodes = pd.DataFrame(adata.obs.index)
barcodes.to_csv('barcodes.tsv', sep='\t', index=False, header=False)