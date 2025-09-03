import scanpy as sc
import pandas as pd
import gseapy as gp
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

adata_path      = "subset_14854X18085_6clu_alpha6_stagate.h5ad"
cluster_key     = "domain"
outdir          = Path("marker_gene_results_14854X18085_6clu_alpha6")
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
adata.write(outdir / "subset_stagate_with_DE.h5ad", compression="gzip")

sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

