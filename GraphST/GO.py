import scanpy as sc
import pandas as pd
import numpy as np
import gseapy as gp
from pathlib import Path
from scipy import sparse

h5ad_path = "cluster1_subclustered_2clu/subset_graphst_with_DE.h5ad"
cluster_key = "sub_leiden"
outdir = Path("cluster1_subclustered_2clu")
go_db = "GO_Biological_Process_2023"

# GSEA
# gsea_clusters = {"4"}
# gsea_perm      = 1000

def select_markers(
        de_df: pd.DataFrame,
        padj: float = 0.05,
        lfc_init: float = 0.8,
        lfc_min: float = 0.3,
        pct_in_th: float = 0.05,
        pct_out_th: float = 0.2,
        min_genes: int = 15,
) -> pd.DataFrame:

    lfc = lfc_init
    while lfc >= lfc_min:
        sel = (de_df["pvals_adj"] < padj) & (de_df["logfoldchanges"] > lfc)
        cand = de_df.loc[sel].copy()
        if len(cand) >= min_genes:
            return cand
        lfc -= 0.2
    cand["pct_in"] = cand["names"].apply(lambda g: pct_in(g))
    cand["pct_out"] = cand["names"].apply(lambda g: pct_out(g))
    sel = (cand["pct_in"] > pct_in_th) & (cand["pct_out"] < pct_out_th)
    return cand.loc[sel]

adata = sc.read_h5ad(h5ad_path)

def pct_in(gene: str) -> float:
    mask = adata.obs[cluster_key] == current_cl
    X = adata[mask, [gene]].X
    if sparse.issparse(X):
        X = X.toarray().ravel()
    return np.mean(X > 0)


def pct_out(gene: str) -> float:
    mask = adata.obs[cluster_key] != current_cl
    X = adata[mask, [gene]].X
    if sparse.issparse(X):
        X = X.toarray().ravel()
    return np.mean(X > 0)


for current_cl in adata.obs[cluster_key].cat.categories:
    print(f"\nCluster {current_cl}")

    de = sc.get.rank_genes_groups_df(adata, group=str(current_cl)).head(3000)

    # if current_cl in gsea_clusters:
    #     cl_dir = outdir / f"GO_cluster_{current_cl}"
    #     cl_dir.mkdir(parents=True, exist_ok=True)
    #
    #     rnk = de.sort_values("logfoldchanges", ascending=False)[["names", "logfoldchanges"]]
    #
    #     pre_res = gp.prerank(
    #         rnk=rnk,
    #         gene_sets=go_db,
    #         organism="Human",
    #         outdir=str(cl_dir),
    #         permutation_num=gsea_perm,
    #         seed=42,
    #         processes=4,
    #     )
    #     pre_res.res2d.to_csv(cl_dir / "GSEA_BP_results.csv", index=False)
    #     continue

    markers = select_markers(de)
    if markers.empty:
        print("No markers found")
        continue

    cl_dir = outdir / f"GO_cluster_{current_cl}"
    cl_dir.mkdir(parents=True, exist_ok=True)
    markers.to_csv(cl_dir / "markers.csv", index=False)
    print(f"saved")

    gene_list = markers["names"].tolist()
    enr = gp.enrichr(
        gene_list=gene_list,
        organism="Human",
        gene_sets=go_db,
        outdir=str(cl_dir),
        cutoff=0.05,
        no_plot=True,
    )
    enr.res2d.to_csv(cl_dir / "GO_BP_enrichr.csv", index=False)
    print(f"saved")

