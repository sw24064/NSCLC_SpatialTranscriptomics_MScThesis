import scanpy as sc
import pandas as pd
import squidpy as sq
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

# --- 输入文件 ---
# adata_path = r"F:\24-25\Final_project\OneDrive_2_2025-5-29\P055_SpatialTranscriptomics\Data\GraphST-main\GraphST-main\subset_14854X18085_leiden_preprocessing_pca_7clu_dim64_radius50_ngbrs15_graphst.h5ad"
adata_path = r"F:\24-25\Final_project\OneDrive_2_2025-5-29\P055_SpatialTranscriptomics\Data\STAGATE-main\STAGATE-main\subset_14854X18085_6clu_alpha6_stagate.h5ad"
adata = sc.read_h5ad(adata_path)

cluster_key = "domain"

custom_annotation = {
    "0": "Epithelial-like",
    "1": "T cells / B cells",
    "2": "Macrophages",
    "3": "Macrophages",
    "4": "Macrophages",
    "5": "Macrophages"
}

adata.obs["celltype_3class"] = adata.obs["domain"].map(custom_annotation)

adata.obs["celltype_3class"] = adata.obs["celltype_3class"].astype("category")

color_map = {
    "Macrophages": "#E64B35",
    "T cells / B cells": "#4DBBD5",
    "Epithelial-like": "#FFD700"
}
adata.uns["celltype_3class_colors"] = [color_map[c] for c in adata.obs["celltype_3class"].cat.categories]

library_id = list(adata.uns["spatial"].keys())[0]
ax = sq.pl.spatial_scatter(
    adata,
    color="celltype_3class",
    library_id=library_id,
    size=1.5,
    img_alpha=0.5,
    legend_loc="right margin",
    return_ax=True
)
ax.set_title("SingleR annotation", fontsize=14)
plt.show()

