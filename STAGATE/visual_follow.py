import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('TkAgg')

adata = sc.read_h5ad("subset_14854X18085_6clu_alpha6_stagate.h5ad")


# # plt.rcParams["figure.figsize"] = (6, 5)
# sc.pl.umap(adata, color='domain', title='STAGATE Clustering (UMAP)')
# sc.pl.umap(adata, color="domain", legend_loc="on data", title="STAGATE Clustering (UMAP)")
#
# # sc.pl.spatial(
# #     adata,
# #     img_key="hires",
# #     color="domain",
# #     title="STAGATE Spatial Clustering",
# #     show=True
# # )
#
# sc.pl.spatial(adata, color="domain", title="Spatial Clustering(STAGATE)", legend_loc="right margin")


# adata.obs["n_counts"] = np.array(adata.X.sum(axis=1)).ravel()
#
# ax = sc.pl.violin(
#     adata,
#     keys='n_counts',
#     groupby='domain',
#     show=False
# )
#
# plt.title("Total UMI counts per Leiden cluster(STAGATE)")
# plt.xlabel("Cluster")
# plt.ylabel("n_counts")
#
# cluster_stats = adata.obs.groupby("domain")["n_counts"].agg(["mean", "median"])
#
# for i, (mean_val, median_val) in enumerate(zip(cluster_stats["mean"], cluster_stats["median"])):
#     plt.hlines(mean_val, i - 0.3, i + 0.3, colors="red", linestyles="--", label="Mean" if i == 0 else "")
#     plt.hlines(median_val, i - 0.3, i + 0.3, colors="blue", linestyles="-", label="Median" if i == 0 else "")
#
# plt.legend()
#
# plt.savefig("violin_n_counts_stagate.png", dpi=600, bbox_inches="tight")
# plt.show()

# annotation
from pandas.api.types import CategoricalDtype
adata = sc.read_h5ad("subset_14854X18085_6clu_alpha6_stagate.h5ad")


anno_map = {
    "0": "Secretory macrophage-like cells",
    "1": "Ig⁺ macrophage-like cells",
    "2": "DC-like secretory cells",
    "3": "Metabolic CAF-like cells",
    "4": "Plasma-like B cell",
    "5": "Activated T cell",
}
adata.obs["cluster_anno"] = adata.obs["domain"].astype(str).map(anno_map)

if hasattr(adata.obs["domain"].dtype, "categories"):
    domain_cats = [str(c) for c in adata.obs["domain"].cat.categories]
else:
    domain_cats = sorted(adata.obs["domain"].astype(str).unique(), key=lambda x: int(x))

cluster_order = [anno_map[c] for c in domain_cats]
adata.obs["cluster_anno"] = adata.obs["cluster_anno"].astype(
    CategoricalDtype(categories=cluster_order, ordered=True)
)

if "domain_colors" in adata.uns:
    dom_color_map = dict(zip(domain_cats, adata.uns["domain_colors"]))
    adata.uns["cluster_anno_colors"] = [dom_color_map[c] for c in domain_cats]

sc.pl.spatial(
    adata,
    img_key="hires",
    color="cluster_anno",
    title="STAGATE annotation",
    show=False
)
plt.savefig("stagate_anno.png", dpi=600, bbox_inches="tight")