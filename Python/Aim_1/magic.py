import magic
import scprep

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import feather

data = scprep.io.load_csv("raw.counts.csv", index_col=0)
data = data.transpose()
data.head()

approx_magic_op = magic.MAGIC(solver="approximate")
data_magic = approx_magic_op.fit_transform(data, genes='all_genes')

data_magic.to_csv("exp_magic.csv")

#feather.write_dataframe(data_magic, "exp_magic.feather")