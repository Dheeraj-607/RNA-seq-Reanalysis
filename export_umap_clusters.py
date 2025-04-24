# export_umap_clusters.py

import scanpy as sc
import os
import json
import numpy as np


adata = sc.read("processed_data_output.h5ad")
output_json = "data/umap_clusters.json"

# Ensure coordinates and cluster
if 'X_umap' not in adata.obsm:
    sc.tl.umap(adata)
if 'leiden' not in adata.obs:
    sc.tl.leiden(adata)

records = []
genes = list(adata.var_names)

# Convert .X to dense if needed (sparse to array)
X = adata.X.toarray() if not isinstance(adata.X, (list, np.ndarray)) else adata.X

# Loop through each cell
for i in range(adata.n_obs):
    cell = {
        "x": float(adata.obsm['X_umap'][i, 0]),
        "y": float(adata.obsm['X_umap'][i, 1]),
        "cell_id": adata.obs_names[i],
        "cluster": str(adata.obs['leiden'][i]),
        "gene_expression": {
            gene: float(X[i, j]) for j, gene in enumerate(genes)
        }
    }
    records.append(cell)

os.makedirs("data", exist_ok=True)
with open(output_json, "w") as f:
    json.dump(records, f, indent=2)

print(f"✅ Exported all genes to: {output_json}")
