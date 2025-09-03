# import scanpy as sc
# import pandas as pd
# import numpy as np
# import gseapy as gp
# from pathlib import Path
# from scipy import sparse
# import os
#
#
# h5ad_path = "marker_gene_results_14854X18085_6clu_alpha6/subset_stagate_with_DE.h5ad"
# cluster_key = "domain"
# outdir = Path("marker_gene_results_14854X18085_6clu_alpha6")
# go_db = "GO_Biological_Process_2023"
#
# # gsea_clusters = {"4"}
# # gsea_perm      = 1000
#
# def select_markers(
#         de_df: pd.DataFrame,
#         padj: float = 0.05,
#         lfc_init: float = 0.8,
#         lfc_min: float = 0.3,
#         pct_in_th: float = 0.05,
#         pct_out_th: float = 0.2,
#         min_genes: int = 15,
# ) -> pd.DataFrame:
#     lfc = lfc_init
#     while lfc >= lfc_min:
#         sel = (de_df["pvals_adj"] < padj) & (de_df["logfoldchanges"] > lfc)
#         cand = de_df.loc[sel].copy()
#         if len(cand) >= min_genes:
#             return cand
#         lfc -= 0.2
#     cand["pct_in"] = cand["names"].apply(lambda g: pct_in(g))
#     cand["pct_out"] = cand["names"].apply(lambda g: pct_out(g))
#     sel = (cand["pct_in"] > pct_in_th) & (cand["pct_out"] < pct_out_th)
#     return cand.loc[sel]
#
# adata = sc.read_h5ad(h5ad_path)
# adata.var_names_make_unique()
#
# def pct_in(gene: str) -> float:
#     mask = adata.obs[cluster_key] == current_cl
#     X = adata[mask, [gene]].X
#     if sparse.issparse(X):
#         X = X.toarray().ravel()
#     return np.mean(X > 0)
#
#
# def pct_out(gene: str) -> float:
#     mask = adata.obs[cluster_key] != current_cl
#     X = adata[mask, [gene]].X
#     if sparse.issparse(X):
#         X = X.toarray().ravel()
#     return np.mean(X > 0)
#
# for current_cl in adata.obs[cluster_key].cat.categories:
#     print(f"\nCluster {current_cl}")
#
#     de = sc.get.rank_genes_groups_df(adata, group=str(current_cl)).head(3000)
#
#     # if current_cl in gsea_clusters:
#     #     cl_dir = outdir / f"GO_cluster_{current_cl}"
#     #     cl_dir.mkdir(parents=True, exist_ok=True)
#     #
#     #     rnk = de.sort_values("logfoldchanges", ascending=False)[["names", "logfoldchanges"]]
#     #     rnk["logfoldchanges"] += np.random.normal(0, 1e-6, size=rnk.shape[0])
#     #
#     #     pre_res = gp.prerank(
#     #         rnk=rnk,
#     #         gene_sets=go_db,
#     #         organism="Human",
#     #         outdir=str(cl_dir),
#     #         permutation_num=gsea_perm,
#     #         seed=42,
#     #         processes=4,
#     #     )
#     #     pre_res.res2d.to_csv(cl_dir / "GSEA_BP_results.csv", index=False)
#     #     print(f"saved")
#     #     continue
#     markers = select_markers(de)
#     if markers.empty:
#         print("skip")
#         continue
#
#     cl_dir = outdir / f"GO_cluster_{current_cl}"
#     cl_dir.mkdir(parents=True, exist_ok=True)
#     markers.to_csv(cl_dir / "markers.csv", index=False)
#     print(f"saved")

#     gene_list = markers["names"].tolist()
#     enr = gp.enrichr(
#         gene_list=gene_list,
#         organism="Human",
#         gene_sets=go_db,
#         outdir=str(cl_dir),
#         cutoff=0.05,
#         no_plot=True,
#     )
#     enr.res2d.to_csv(cl_dir / "GO_BP_enrichr.csv", index=False)
#     print(f"saved")





# bar chart
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('TkAgg')

file_path = "marker_gene_results_14854X18085_6clu_alpha6/GO_cluster_5/GO_BP_enrichr.csv"
cluster_id = "5"

df = pd.read_csv(file_path)
df = df[df["Adjusted P-value"] < 0.05]
df = df[df["Combined Score"] > 30]

df = df.sort_values("Adjusted P-value", ascending=True).head(20)
df = df.sort_values("Combined Score", ascending=False).head(10)



plt.figure(figsize=(8, 6))
sns.barplot(data=df, x="Combined Score", y="Term", palette="magma")

plt.title(f"Top 10 GO Terms for Cluster {cluster_id} (STAGATE)")
plt.xlabel("Combined Score")
plt.ylabel("GO Term")
plt.tight_layout()
plt.savefig("Top10_GO_Cluster5_stagate.png", dpi=600, bbox_inches="tight")
plt.show()


