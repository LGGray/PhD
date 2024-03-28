import json
import zlib
import base64
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp

# collect SCENIC AUCell output
lf = lp.connect( "SCENIC/pbmc.female.SCENIC.loom", mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
pd.DataFrame.to_csv(auc_mtx, "SCENIC/auc_mtx.csv")
lf.close()