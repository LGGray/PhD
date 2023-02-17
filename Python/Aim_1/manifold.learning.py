import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

# Load the single-cell RNAseq data
data = pd.read_csv("data.csv")

# Extract the gene expression values from the data
X = np.array(data.iloc[:, 1:])

# Compute the t-SNE embedding of the data
tsne = TSNE(n_components=2, perplexity=30, random_state=0)
X_tsne = tsne.fit_transform(X)

# Extract the loadings of each gene
gene_loadings = tsne.embedding_.T

# Identify the top genes that are driving the separation of the clusters
n_top_genes = 20
top_genes_idx = np.argsort(gene_loadings[0])[::-1][:n_top_genes]
top_genes = data.columns[1:][top_genes_idx]

# Plot the t-SNE embedding with the top genes highlighted
plt.scatter(X_tsne[:, 0], X_tsne[:, 1], s=5)
plt.scatter(X_tsne[:, 0][top_genes_idx], X_tsne[:, 1][top_genes_idx], s=50, c='r')
plt.title("t-SNE visualization of single-cell RNAseq data with top genes highlighted")
plt.show()
