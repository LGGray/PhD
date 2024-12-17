import scanpy as sc
import pandas as pd
from scipy.sparse import save_npz

# Read in .h5ad
adata = sc.read_h5ad("4532eea4-24b7-461a-93f5-fe437ee96f0a.h5ad")

# Subset for female sex and managed or na disease state
adata_female = adata[(adata.obs['sex'] == 'female') & ((adata.obs['disease_state'] == 'managed') | (adata.obs['disease_state'] == 'na')), :]

# Export sparse matrix
save_npz('exp_counts_female_managed.npz', adata_female.X)

# Export metadata
metadata = pd.DataFrame(adata_female.obs)
metadata.to_csv('metadata.tsv', sep='\t', index=True)

# Export feature
features = pd.DataFrame(adata_female.var)
features['gene_id'] = features.index
features = features[['gene_id', 'feature_name', 'feature_biotype']]
features.to_csv('features.tsv', sep='\t', index=False, header=False)

# Export barcodes
barcodes = pd.DataFrame(adata_female.obs.index)
barcodes.to_csv('barcodes.tsv', sep='\t', index=False, header=False)

#######

# Subset for female sex and flare or na disease state
adata_flare = adata[(adata.obs['sex'] == 'female') & ((adata.obs['disease_state'] == 'flare') | (adata.obs['disease_state'] == 'na')), :]

# Export sparse matrix
save_npz('FLARE/exp_counts_female_flare.npz', adata_flare.X)

# Export metadata
metadata = pd.DataFrame(adata_flare.obs)
metadata.to_csv('FLARE/metadata.tsv', sep='\t', index=True)

# Export feature
features = pd.DataFrame(adata_flare.var)
features['gene_id'] = features.index
features = features[['gene_id', 'feature_name', 'feature_biotype']]
features.to_csv('FLARE/features.tsv', sep='\t', index=False, header=False)

# Export barcodes
barcodes = pd.DataFrame(adata_flare.obs.index)
barcodes.to_csv('FLARE/barcodes.tsv', sep='\t', index=False, header=False)
