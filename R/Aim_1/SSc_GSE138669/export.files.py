import scanpy as sc
from scipy.sparse import save_npz
import pandas as pd
import glob
import anndata



# Get the list of .h5 files in the current directory
files = glob.glob("h5.files/*.h5")
samples = [file.split('/')[1].split('_')[1].split('raw')[0] for file in files]
# Create an empty list to store the AnnData objects
adata_list = []
# Loop through the files and read them using scanpy.read_10x_h5
for file in files:
    print(file)
    adata = sc.read_10x_h5(file)
    # Add the sample name to the AnnData object
    adata.obs['individual'] = samples[files.index(file)]
    # Add the gene symbol to the variable metadata
    adata.var['gene_symbol'] = adata.var.index
    # Add ENSG ID to the variable index
    adata.var.index = adata.var['gene_ids']
    # Append the AnnData object to the list
    adata_list.append(adata)
# Concatenate the AnnData objects in the list
adata = anndata.concat(adata_list, index_unique='-')

# Export sparse matrix
save_npz('exp_counts.npz', adata.X)
metadata = pd.DataFrame(adata.obs)
metadata['batch'] = metadata['individual']
metadata.to_csv('cell_batch.tsv', sep='\t', index=True)
# Align the variable metadata
features = pd.DataFrame(adata.var.align(adata_list[0].var, join='inner', axis=0)[1])
features = features[['gene_ids', 'gene_symbol', 'feature_types']]
features.to_csv('features.tsv', sep='\t', index=False, header=False)
barcodes = pd.DataFrame(adata.obs.index)
barcodes.to_csv('barcodes.tsv', sep='\t', index=False, header=False)