import scvelo as scv
import matplotlib.pyplot as plt

with open('../processed/merged_adata.pkl', 'rb') as f:
        merged_adata = pickle.load(f)

scv.set_figure_params()

scv.pl.velocity_embedding_stream(merged_adata, basis='umap', color='seurat_clusters', figsize=(10, 10), save='../figs/paper/scv1.svg')

with open('../processed/microglial_adata.pkl', 'rb') as f:
        microglial_adata = pickle.load(f)

scv.pl.velocity_embedding_stream(microglial_adata, basis='umap', color='seurat_clusters', figsize=(10, 10), save='../figs/paper/scv2.svg')

scv.tl.velocity_pseudotime(microglial_adata)
scv.pl.scatter(microglial_adata, color = "velocity_pseudotime", color_map = "viridis", show = False, figsize = (10, 11), save='../figs/paper/scv3.svg')
