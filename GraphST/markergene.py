import scanpy as sc
import pandas as pd
import gseapy as gp
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

adata_path      = "cluster1_subclustered_2clu.h5ad"
cluster_key     = "sub_leiden"
outdir          = Path("cluster1_subclustered_2clu")
outdir.mkdir(exist_ok=True)

adata = sc.read_h5ad(adata_path)

sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon")

all_marker = sc.get.rank_genes_groups_df(adata, None)
all_marker.to_csv(outdir / "rank_genes_groups_all_clusters.csv", index=False)

clusters = adata.obs[cluster_key].astype(str).unique()
for cl in clusters:
    cluster_dir = outdir / f"GO_cluster_{cl}"
    cluster_dir.mkdir(parents=True, exist_ok=True)

    df = sc.get.rank_genes_groups_df(adata, group=cl)
    gene_list = df["names"].tolist()

    enr = gp.enrichr(
        gene_list=gene_list,
        organism="Human",
        gene_sets="GO_Biological_Process_2023",
        outdir=str(cluster_dir),
        cutoff=0.05,
        no_plot=True,
    )
    enr.res2d.to_csv(cluster_dir / "GO_BP_enrichr.csv", index=False)

adata.raw = adata
expr_df = adata.raw.to_adata().to_df()
expr_df.to_csv(outdir / "expression_matrix.csv")
adata.obs[[cluster_key]].to_csv(outdir / "cluster_labels.csv")
adata.write(outdir / "subset_graphst_with_DE.h5ad", compression="gzip")




# vocano

# import pandas as pd
# import numpy as np
# import scanpy as sc
# from adjustText import adjust_text
# import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.use('TkAgg')
#
# h5ad_path  = r"F:\24-25\Final_project\OneDrive_2_2025-5-29\P055_SpatialTranscriptomics\Data\GraphST-main\GraphST-main\marker_gene_results_14854X18085_7clu_15ngbrs\subset_graphst_with_DE.h5ad"
# marker_csv = r"F:\24-25\Final_project\OneDrive_2_2025-5-29\P055_SpatialTranscriptomics\Data\GraphST-main\GraphST-main\marker_gene_results_14854X18085_7clu_15ngbrs\GO_cluster_0\markers.csv"
# cluster_key, target_grp = "leiden", "0"
#
# adata = sc.read_h5ad(h5ad_path)
# rgg = adata.uns["rank_genes_groups"]
# groups = list(rgg["names"].dtype.names)
# assert target_grp in groups, f"{target_grp=} is not in {groups}"
#
# def col_or_nan(key):
#     return rgg[key][target_grp] if key in rgg else np.full(len(rgg["names"][target_grp]), np.nan)
#
# de = pd.DataFrame({
#     "gene":          rgg["names"][target_grp].astype(str),
#     "logFC":         col_or_nan("logfoldchanges"),
#     "padj":          col_or_nan("pvals_adj"),
#     "pct_in_group":  col_or_nan("pts"),
#     "pct_out_group": col_or_nan("pts_rest"),
# }).dropna(subset=["logFC","padj"]).set_index("gene")
#
# mk = pd.read_csv(marker_csv)
# name_col = [c for c in mk.columns if c.lower() in {"gene","genes","names","name"}]
# assert len(name_col)>=1, f"marker.csv cannot find：{mk.columns}"
# mk_genes = mk[name_col[0]].astype(str).unique()
#
# plot_df = de.loc[de.index.intersection(mk_genes)].copy()
# plot_df["neglog10padj"] = -np.log10(plot_df["padj"].clip(lower=1e-300))
#
# top10 = plot_df.nsmallest(10, "padj")
#
# plt.figure(figsize=(8.5,7.5))
#
# others = de.index.difference(plot_df.index)
# if len(others):
#     plt.scatter(de.loc[others, "logFC"], -np.log10(de.loc[others,"padj"].clip(lower=1e-300)),
#                 s=10, c="#e6e6e6", alpha=0.5, linewidths=0)
#
# size_base = 600
# sizes_all = (plot_df["pct_in_group"].fillna(0.15) * size_base).clip(60, 900)
# plt.scatter(plot_df["logFC"], plot_df["neglog10padj"],
#             s=30, c="#ff944d", alpha=0.5, linewidths=0, label="marker (all)")
#
# sizes_top = (top10["pct_in_group"].fillna(0.15) * (size_base*1.2)).clip(80, 1100)
# plt.scatter(top10["logFC"], top10["neglog10padj"],
#             s=sizes_top, c="#d62728", alpha=0.95, linewidths=0, label="Top10 (by FDR)")
#
# padj_thr, fc_thr = 0.05, 1.0
# plt.axhline(-np.log10(padj_thr), ls="--", c="black", lw=1)
# plt.axvline(fc_thr,             ls="--", c="black", lw=1)
#
# texts = []
# for gene, row in top10.iterrows():
#     texts.append(
#         plt.text(
#             row.logFC,
#             row.neglog10padj,
#             gene,
#             fontsize=12,
#             ha="center",
#             va="bottom"
#         )
#     )
#
# adjust_text(
#     texts,
#     arrowprops=dict(arrowstyle="-", color="black", lw=0.5),
#     expand_points=(1.2, 1.6),
#     expand_text=(1.2, 1.6),
#     force_points=0.3,
#     force_text=0.5,
#     only_move={'points':'y', 'text':'xy'}
# )
#
# plt.title(f"Cluster {target_grp} (N={len(plot_df)})", fontsize=18)
# plt.xlabel("log2 Fold Change (cluster vs rest)", fontsize=14)
# plt.ylabel("-log10(FDR)", fontsize=14)
# plt.legend(loc="upper left", frameon=False, fontsize=12)
# plt.tight_layout()
# plt.show()