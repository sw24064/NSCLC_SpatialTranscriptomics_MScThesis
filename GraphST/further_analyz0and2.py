import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

adata = sc.read_h5ad("subset_14854X18085_leiden_preprocessing_pca_7clu_dim64_radius50_ngbrs15_graphst.h5ad")

cluster_col = "leiden"
target_cluster = "1"

adata.obs["highlight_cluster"] = adata.obs["leiden"].astype(str)
adata.obs["highlight_cluster"] = adata.obs["highlight_cluster"].apply(
    lambda x: x if x == target_cluster else "other"
)

# sc.pl.spatial(
#     adata,
#     color="highlight_cluster",
#     img_key="hires",
#     palette=["crimson", "lightgray"],
#     size=1.3,
#     title=f"Spatial distribution of cluster {target_cluster}",
#     show=True
# )

adata_sub = adata[adata.obs["leiden"] == "1", :].copy()

sc.pp.pca(adata_sub)
sc.pp.neighbors(adata_sub, n_neighbors=5)
for res in [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
    sc.tl.leiden(adata_sub, resolution=res, key_added=f"leiden_{res}")
    n = adata_sub.obs[f"leiden_{res}"].nunique()
    print(f"resolution={res:.2f} → {n} clusters")
sc.tl.leiden(adata_sub, resolution=0.4)
# adata_sub.write("cluster0_subclustered_2clu.h5ad")

sc.tl.umap(adata_sub)
sc.pl.umap(
    adata_sub,
    color="leiden",
    title="UMAP of subclusters (from cluster 1)",
    size=30
)

sc.pl.spatial(
    adata_sub,
    color="leiden",
    img_key="hires",
    size=1.3,
    title="Spatial distribution of subclusters (from cluster 1)"
)



# sc.tl.rank_genes_groups(
#     adata_sub,
#     groupby="leiden",
#     method="wilcoxon",
#     reference="rest"
# )
#
#
# sc.pl.rank_genes_groups_dotplot(
#     adata_sub,
#     n_genes=10,
#     groupby="leiden",
#     standard_scale="var",
#     dendrogram=True
# )
#




