import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from pandas.api.types import CategoricalDtype
import matplotlib.patches as mpatches
from sklearn.metrics import silhouette_score
import scipy.sparse

import matplotlib
matplotlib.use('TkAgg')


adata = sc.read_h5ad("subset_14854X18085_leiden_preprocessing_pca_7clu_dim64_radius50_ngbrs15_graphst.h5ad")


# # violin
# adata.obs["n_counts"] = np.array(adata.X.sum(axis=1)).ravel()
# ax = sc.pl.violin(
#     adata,
#     keys='n_counts',
#     groupby='leiden',
#     show=False
# )
#
# plt.title("Total UMI counts per Leiden cluster")
# plt.xlabel("Cluster")
# plt.ylabel("n_counts")
#
# cluster_stats = adata.obs.groupby("leiden")["n_counts"].agg(["mean", "median"])
#
# for i, (mean_val, median_val) in enumerate(zip(cluster_stats["mean"], cluster_stats["median"])):
#     plt.hlines(mean_val, i - 0.3, i + 0.3, colors="red", linestyles="--", label="Mean" if i == 0 else "")
#     plt.hlines(median_val, i - 0.3, i + 0.3, colors="blue", linestyles="-", label="Median" if i == 0 else "")
#
# plt.legend()
#
# plt.savefig("violin_n_counts.png", dpi=600, bbox_inches="tight")
# plt.show()
# print(adata)
# print("emb:", adata.obsm['emb'].shape)
# print("emb_pca:", adata.obsm['emb_pca'].shape)
#
# print(adata.uns_keys())
# print(adata.obs["leiden"].value_counts())
#
# print(adata.obs.columns)



# spatial
fig = sc.pl.spatial(
    adata,
    img_key="hires",
    color="leiden",
    title="GraphST Spatial Clustering",
    show=False,
    return_fig=True
)
plt.savefig("graphst_SPATIAL.png", dpi=600, bbox_inches="tight")
# UMAP
sc.pp.neighbors(adata, use_rep='emb_pca', n_neighbors=10)
sc.tl.umap(adata)
sc.pl.umap(adata, color='leiden', title='GraphST Clustering (UMAP)')



# Final annotation
sub0 = sc.read_h5ad("cluster0_subclustered_2clu.h5ad")
sub1  = sc.read_h5ad("cluster1_subclustered_2clu.h5ad")

sub0col = "sub_leiden" if "sub_leiden" in sub0.obs.columns else "leiden"
sub1col = "sub_leiden" if "sub_leiden" in sub1.obs.columns else "leiden"
sub0.obsm["spatial"] = adata.obsm["spatial"][[adata.obs_names.get_loc(x) for x in sub0.obs_names], :]
sub1.obsm["spatial"] = adata.obsm["spatial"][[adata.obs_names.get_loc(x) for x in sub1.obs_names], :]

adata.obs["cluster_anno"] = adata.obs["leiden"].astype(str).replace({"3": "1"})
anno_map = {
    "0": "CAFs",
    "1": "TLS-associated B cells / Plasma cells",
    "2": "TAMs",
    "4": "Proliferating TAMs with CAF/ECM signals",
    "5": "Inflammatory CAF–TAMs",
    "6": "SPP1⁺ TAMs"
}
adata.obs["cluster_anno"] = adata.obs["cluster_anno"].map(anno_map)

cluster_order = [
    "CAFs",
    "TLS-associated B cells / Plasma cells",
    "TAMs",
    "Proliferating TAMs with CAF/ECM signals",
    "Inflammatory CAF–TAMs",
    "SPP1⁺ TAMs",
]
major_colors = [
    "#377eb8",
    "#ff7f00",
    "#4daf4a",
    "#984ea3",
    "#a65628",
    "#f781bf",
]

adata.obs["cluster_anno"] = adata.obs["cluster_anno"].astype(
    CategoricalDtype(categories=cluster_order, ordered=True)
)
adata.uns["cluster_anno_colors"] = major_colors

adata.obs.loc[sub0.obs_names, "CAF_sub"] = sub0.obs[sub0col]
adata.obs["CAF_sub"] = adata.obs["CAF_sub"].map({"0": "iCAFs", "1": "ECM-remodeling CAFs"})
caf_order = ["iCAFs", "ECM-remodeling CAFs"]
adata.obs["CAF_sub"] = adata.obs["CAF_sub"].astype(CategoricalDtype(categories=caf_order, ordered=True))
caf_colors = ["#99ccff", "#1f4e79"]

adata.obs.loc[sub1.obs_names, "TLS_sub"] = sub1.obs[sub1col]
adata.obs["TLS_sub"] = adata.obs["TLS_sub"].map({"0": "Plasma-like B cells", "1": "T cell–dominant immune niche with B cell"})
tls_order  = ["Plasma-like B cells", "T cell–dominant immune niche with B cell"]
tls_colors = ["#8dd3c7", "#ffffb3"]
adata.obs["TLS_sub"] = adata.obs["TLS_sub"].astype(
    CategoricalDtype(categories=tls_order, ordered=True)
)


fig = sc.pl.spatial(
    adata,
    img_key="hires",
    color="cluster_anno",
    title="GraphST annotation",
    show=False,
    return_fig=True
)
ax = fig.axes[0]

caf_mask = adata.obs["cluster_anno"] == "CAFs"
sc.pl.spatial(
    adata[caf_mask],
    img_key="hires",
    color="CAF_sub",
    palette=caf_colors,
    ax=ax,
    show=False
)

tls_mask = adata.obs["cluster_anno"] == "TLS-associated B cells / Plasma cells"
sc.pl.spatial(
    adata[tls_mask],
    img_key="hires",
    color="TLS_sub",
    palette=tls_colors,
    ax=ax,
    show=False
)

# ax.invert_yaxis()
ax.set_title("Annotation")

main_handles = [mpatches.Patch(color=c, label=l) for c, l in zip(major_colors, cluster_order)]
caf_handles  = [mpatches.Patch(color=c, label=l) for c, l in zip(caf_colors, caf_order)]
tls_handles  = [mpatches.Patch(color=c, label=l) for c, l in zip(tls_colors, tls_order)]
leg1 = ax.legend(handles=main_handles, bbox_to_anchor=(1.05, 1.0), loc="upper left", title="Major clusters")
leg2 = ax.legend(handles=caf_handles,  bbox_to_anchor=(1.05, 0.55), loc="upper left", title="CAFs subclusters")
leg3 = ax.legend(handles=tls_handles,  bbox_to_anchor=(1.05, 0.25), loc="upper left", title="TLS-associated B cells / Plasma cells subclusters")
ax.add_artist(leg1); ax.add_artist(leg2)

plt.savefig("final_annotation_sub.png", dpi=600, bbox_inches="tight")
plt.show()
