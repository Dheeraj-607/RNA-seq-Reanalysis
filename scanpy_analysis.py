# scanpy_analysis.py

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os
import gzip

# Config
data_dir = r"C:/Users/Dell\Desktop/steered-project_Dheeraj/GSE67835_RAW"
output_csv = "gene_counts_matrix_output.csv"
output_h5ad = "processed_data_output.h5ad"
umap_png = "figures/umap_clusters_member2.png"

# Load .csv.gz expression files
files = [f for f in os.listdir(data_dir) if f.endswith(".csv.gz")]
all_cells = []

for f in files:
    file_path = os.path.join(data_dir, f)
    try:
        with gzip.open(file_path, 'rt') as fh:
            df = pd.read_csv(fh, index_col=0, header=None, sep='\t')
        if df.shape[0] == 0 or df.shape[1] != 1:
            print(f"Skipping malformed: {f}")
            continue
        df.columns = [f.split("_")[0]]
        all_cells.append(df)
    except Exception as e:
        print(f"Error reading {f}: {e}")

if not all_cells:
    raise RuntimeError("No valid input files were loaded.")

combined_df = pd.concat(all_cells, axis=1)
combined_df.to_csv(output_csv)

# Create AnnData and preprocess
adata = sc.AnnData(combined_df.T)
adata.var.index.name = None
adata.var_names_make_unique()

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# Save output
sc.pl.umap(adata, color='leiden', save="_clusters_image.png")
adata.write(output_h5ad)

print("✅ Analysis complete.")
