import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import cm
import matplotlib.colors as mcolors
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import random
from scipy import sparse
import matplotlib
matplotlib.use('TkAgg')

'''raw'''
# adata = sc.read_h5ad(
#     r"F:\24-25\Final_project\OneDrive_2_2025-5-29\P055_SpatialTranscriptomics\Data\GraphST-main\GraphST-main\GraphST_trained_withSig.h5ad"
# )
#
# adata_full = sc.read_h5ad(
#     r"F:\24-25\Final_project\OneDrive_2_2025-5-29\P055_SpatialTranscriptomics\Data\SpatialProject\GraphST-main\subset_14854X18085.h5ad"
# )
#
# adata_full = adata_full[adata.obs_names, :]
#
# adata.raw = adata_full
#
# out = r"F:\24-25\Final_project\GraphST_withRaw.h5ad"
# adata.write(out)


adata = sc.read_h5ad(r"F:\24-25\Final_project\GraphST_withRaw.h5ad")
sig_df = pd.read_csv(r"C:\Users\FREYA\Desktop\signatures.csv", header=0)
sig_dict = {c: sig_df[c].dropna().str.upper().tolist() for c in sig_df.columns}


sig_name   = "CAF"
score_col  = f"sig_{sig_name}"
gene_list  = [g.strip().upper() for g in sig_dict[sig_name] if isinstance(g, str) and g.strip() != '']
print(f"{sig_name} gene list:", len(gene_list))
print("gene list:", gene_list[:10])

'''heatmap'''
# cluster_col = "domain"
# adata_raw = adata.raw.to_adata()
# adata_raw.var_names_make_unique()
# adata_raw.var_names = adata_raw.var_names.str.upper()
#
# gene_list_use = [g for g in gene_list if g in adata_raw.var_names]

# if len(gene_list_use) == 0:
#     raise ValueError("no gene")
#
# X = adata_raw[:, gene_list_use].X
# X = X.toarray() if sparse.issparse(X) else X
#
# X_df = pd.DataFrame(X, columns=gene_list_use, index=adata.obs_names)
# X_df[cluster_col] = adata.obs[cluster_col].values
#
# mean_df = X_df.groupby(cluster_col).mean()
#
# plt.figure(figsize=(0.6 * len(gene_list_use), 4))
# sns.heatmap(mean_df, cmap='viridis')
# plt.title(f"Mean expression per cluster: {sig_name} (GraphST)")
# plt.ylabel("Cluster")
# plt.xlabel("Gene")
# plt.tight_layout()
# plt.savefig("IFNγ_iCAF_graphst.png", dpi=600, bbox_inches="tight")
# plt.show()

'''metric'''

# adata.var_names_make_unique()
# adata.var_names = adata.var_names.str.upper()
#
# var_names_full = adata.raw.var_names if adata.raw is not None else adata.var_names
#
# hits = sum(g in var_names_full for g in gene_list)
# print("accuracy:", hits / len(gene_list))
#
# print("in .var  match:", sum(g in adata.var_names for g in gene_list))
# if adata.raw is not None:
#     print("in .raw  match:", sum(g in adata.raw.var_names for g in gene_list))
# else:
#     print("no raw")

# marker = "MS4A1"
# sc.pl.scatter(adata, x=score_col, y=marker)

# lib_id = list(adata.uns["spatial"].keys())[0]
# vmax99 = np.percentile(adata.obs[score_col], 99)
# sq.pl.spatial_scatter(
#     adata,
#     library_id=lib_id,
#     color=score_col,
#     cmap="magma",
#     img_alpha=0.4,
#     size=1.2,
#     vmax=vmax99,
#     title=sig_name
# )
# plt.show()

'''for single gene'''
lib_id = list(adata.uns["spatial"].keys())[0]

fig, ax = plt.subplots(figsize=(6, 6))

vmax99 = np.percentile(adata.obs[score_col], 99)
sq.pl.spatial_scatter(
    adata,
    library_id=lib_id,
    color=score_col,
    cmap="magma",
    img_alpha=0.4,
    size=1.2,
    vmax=vmax99,
    ax=ax
)

marker_colors = {"CXCL12": "royalblue"}
for gene, col in marker_colors.items():
    vals = adata[:, gene].X
    if hasattr(vals, "toarray"):
        vals = vals.toarray().flatten()
    else:
        vals = vals.flatten()

    vmax_g = np.percentile(vals, 99)
    mask = vals > vmax_g * 0.7

    adata.obs[f"{gene}_highlight"] = np.where(mask, vals, np.nan)

    sq.pl.spatial_scatter(
        adata,
        library_id=lib_id,
        color=f"{gene}_highlight",
        size=1.2,
        cmap=mcolors.ListedColormap([col]),
        ax=ax,
        colorbar=False
    )
ax.set_title("CAF: CXCL12")
patches = [mpatches.Patch(color=col, label=gene) for gene, col in marker_colors.items()]
plt.legend(handles=patches, title="Markers", loc="upper right")
plt.savefig("CAF_CXCL12_1gene.png", dpi=600, bbox_inches="tight")
plt.show()

