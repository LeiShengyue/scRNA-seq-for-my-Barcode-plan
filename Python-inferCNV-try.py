import scanpy as sc
import pandas as pd
import numpy as np
from infercnvpy import infercnv

# Load the data
data1 = sc.read("/path/to/MPS1ip6.h5ad")  # Assuming .h5ad format for Scanpy
data2 = sc.read("/path/to/MPS1ip1.h5ad")

# Preprocess the data
sc.pp.filter_cells(data1, min_genes=200)
sc.pp.filter_cells(data2, min_genes=200)
sc.pp.filter_genes(data1, min_cells=3)
sc.pp.filter_genes(data2, min_cells=3)

# Normalize the data
sc.pp.normalize_total(data1, target_sum=1e4)
sc.pp.normalize_total(data2, target_sum=1e4)
sc.pp.log1p(data1)
sc.pp.log1p(data2)

# Merge the datasets
adata = data1.concatenate(data2, batch_key='batch')

# Find highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# PCA analysis
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# Neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# UMAP visualization
sc.tl.umap(adata)
sc.pl.umap(adata, color=['batch'])

# Differential expression analysis
sc.tl.rank_genes_groups(adata, groupby='batch', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

# Save the processed data
adata.write("/path/to/processed_adata.h5ad")

# -------- CNV Analysis --------

# Define the reference (normal) cells
# Assuming you have annotated your cells with a 'cell_type' column where normal cells are labeled as 'Normal'
normal_cells = adata[adata.obs['cell_type'] == 'Normal']

# Create an InferCNV object using infercnvpy
cnv_obj = infercnv.create_infercnv_object(
    adata=adata,
    reference_key='cell_type',
    reference_cat='Normal',  # Replace 'Normal' with your specific reference label
    window_size=100,  # You can adjust this based on your data
    max_centered_threshold=3,
)

# Run CNV analysis
infercnv.tl.run(cnv_obj)

# Plot the CNV results
infercnv.pl.chromosomal_heatmap(cnv_obj, save="CNV_heatmap.png")

# Optionally save the CNV data for further analysis
cnv_obj.write("/path/to/cnv_results.h5ad")

# View CNV results in a UMAP plot colored by inferred CNV state
sc.pl.umap(adata, color=['cnv_states'])  # 'cnv_states' would be the column added by infercnvpy

# Save the final AnnData object with CNV information
adata.write("/path/to/processed_adata_with_cnv.h5ad")
