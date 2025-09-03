import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

#ROI
x_min, x_max = 8500, 9600
y_min, y_max = 6600, 7999

roi = adata[
    (adata.obsm["spatial"][:, 0] >= x_min) & (adata.obsm["spatial"][:, 0] <= x_max) &
    (adata.obsm["spatial"][:, 1] >= y_min) & (adata.obsm["spatial"][:, 1] <= y_max)
].copy()

print(f"ROI size: {roi.n_obs} spots")

sc.pp.neighbors(roi, use_rep="emb", n_neighbors=8)

for res in [0.1, 0.2, 0.3, 0.4, 0.5]:
    sc.tl.leiden(roi, resolution=res, key_added=f"roi_leiden_{res}")
    n = roi.obs[f"roi_leiden_{res}"].nunique()
    print(f"resolution={res:.2f} → {n} clusters")

sc.tl.leiden(roi, resolution=0.5, key_added="roi_leiden")

sc.tl.umap(roi)
sc.pl.umap(roi, color="roi_leiden", title="ROI clustering (based on GraphST emb)")
sc.pl.spatial(
    roi,
    img_key="hires",
    color="roi_leiden",
    title="GraphST Spatial Clustering",
    show=True
)

