# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import glob

# set working directory
os.chdir('SCENIC')

# Path to loom file
f_loom_path_scenic = "pbmc.loom"

# path to anndata object, which will be updated to store Scanpy results as they are generated below
f_anndata_path = "anndata.h5ad"

# path to pyscenic output
f_pyscenic_output = "SCENIC/pyscenic_output.loom"

# loom output, generated from a combination of Scanpy and pySCENIC results:
f_final_loom = 'SCENIC/pbmc10k_scenic_integrated-output.loom'

# Set plotting parameters
sc.set_figure_params(dpi=150, fontsize=10, dpi_save=600)

### STEP 1: Gene regulatory network inference and generation of co-expression modules with pySCENIC ###
# transcription factors list
f_tfs = "/directflow/SCCGGroupShare/projects/lacgra/SCENIC/allTFs_hg38.txt"

# pyscenic grn {f_loom_path_scenic} {f_tfs} -o adj.csv --num_workers 20

# Read in adjacencies matrix
adjacencies = pd.read_csv("adj.tsv", index_col=False, sep='\t')

### STEP 2-3: Regulon prediction aka cisTarget from CLI ###
f_db_glob = "/directflow/SCCGGroupShare/projects/lacgra/SCENIC/cisTarget_databases/*feather"
f_db_names = ''.join (glob.glob(f_db_glob))

# Motif databases
f_motif_path = 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'

# !pyscenic ctx adj.tsv \
#     {f_db_names} \
#     --annotations_fname {f_motif_path} \
#     --expression_mtx_fname {f_loom_path_scenic} \
#     --output reg.csv \
#     --mask_dropouts \
#     --num_workers 20

## Step 4: Cellular enrichment (AUCell) from CLI ###
# !pyscenic aucell \
#     {f_loom_path_scenic} \
#     reg.csv \
#     --output {f_pyscenic_output} \
#     --num_workers 20

# Visualisation of SCENIC AUC matrix
lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()


