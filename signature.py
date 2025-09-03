import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
import matplotlib.pyplot as plt
import random
from scipy import sparse
import matplotlib
matplotlib.use('TkAgg')
import gseapy as gp


plt.rcParams.update({'font.size': 7})
# RAW_H5AD        = r"F:\24-25\Final_project\OneDrive_2_2025-5-29\P055_SpatialTranscriptomics\Data\SpatialProject\GraphST-main\subset_14854X18085.h5ad"
# H5AD    = r"F:\24-25\Final_project\OneDrive_2_2025-5-29\P055_SpatialTranscriptomics\Data\STAGATE-main\STAGATE-main\subset_14854X18085_6clu_alpha6_stagate.h5ad"
# SIGNATURE_CSV   = r"C:\Users\FREYA\Desktop\signatures.csv"
# OUT_H5AD        = r"F:\24-25\Final_project\OneDrive_2_2025-5-29\P055_SpatialTranscriptomics\Data\STAGATE-main\STAGATE-main\stagate_trained_withSig.h5ad"
#
# adata_full = sc.read_h5ad(RAW_H5AD)
# sc.pp.normalize_total(adata_full, target_sum=1e4)
# sc.pp.log1p(adata_full)
#
# adata_full.var_names_make_unique()
# adata_full.var_names = adata_full.var_names.str.upper()
#
# sig_df = pd.read_csv(SIGNATURE_CSV, header=0)
# sig_dict = {
#     col: sig_df[col].dropna().str.upper().tolist()
#     for col in sig_df.columns
# }
#
# adata_graphst = sc.read_h5ad(H5AD)
#
# if not all(adata_full.obs_names == adata_graphst.obs_names):
#     adata_graphst = adata_graphst[adata_full.obs_names, :]
#
# for name, glist in sig_dict.items():
#     genes_use = [g for g in glist if g in adata_full.var_names]
#     if len(genes_use) < 3:
#         print(f"{name}: match {len(genes_use)} genes")
#         continue
#
#     X = adata_full[:, genes_use].X
#     X = X.toarray() if sparse.issparse(X) else X
#
#     med = np.median(X, axis=0)
#     mad = np.median(np.abs(X - med), axis=0)
#     mad[mad == 0] = 1e-12
#
#     adata_graphst.obs[f"sig_{name}"] = ((X - med) / mad).mean(axis=1)
#

# adata_graphst.write(OUT_H5AD)
# print("Saved:", OUT_H5AD)
#
# # if "X_umap" not in adata_graphst.obsm:
# #     sc.pp.neighbors(adata_graphst, use_rep="emb_pca")
# #     sc.tl.umap(adata_graphst)
# #
# # sig_cols = [f"sig_{n}" for n in sig_dict.keys() if f"sig_{n}" in adata_graphst.obs]
# #
# # sc.pl.umap(
# #     adata_graphst,
# #     color=sig_cols,
# #     cmap="viridis",
# #     ncols=3,
# #     wspace=0.35
# # )
#
# for n, glist in sig_dict.items():
#     missing = [g for g in glist if g not in adata_full.var_names]
#     print(f"{n}: match {len(glist)-len(missing)}/{len(glist)}  miss: {missing[:5]}")






'''plot the mapped-back spatial map'''
adata = sc.read_h5ad(
    r"F:\24-25\Final_project\OneDrive_2_2025-5-29\P055_SpatialTranscriptomics\Data\GraphST-main\GraphST-main\graphst_trained_withSig.h5ad"
)
#
# sig_col = "sig_CAF"
# display_name = sig_col.replace("sig_", "", 1)
#
# LIB = list(adata_graphst.uns["spatial"].keys())[0]
#
# vmax_99 = np.percentile(adata_graphst.obs[sig_col], 99)
#
# sq.pl.spatial_scatter(
#     adata_graphst,
#     library_id=LIB,
#     color=sig_col,
#     cmap="magma",
#     img_alpha=0.4,
#     size=1.2,
#     vmax=vmax_99,
#     frameon=False,
#     title=display_name
# )
# plt.show()


'''check whether the signature meets the criteria'''
sig_df = pd.read_csv(r"C:\Users\FREYA\Desktop\signatures.csv", header=0)
sig_dict = {c: sig_df[c].dropna().str.upper().tolist() for c in sig_df.columns}

sig_name   = "TLS_mature"
score_col  = f"sig_{sig_name}"
gene_list  = sig_dict[sig_name]

adata.var_names_make_unique()
adata.var_names = adata.var_names.str.upper()
hits = sum(g in adata.var_names for g in gene_list)
print("命中率:", hits / len(gene_list))


lm = smf.ols(f'{score_col} ~ C(domain)', data=adata.obs).fit()
anova = anova_lm(lm, typ=2)
eta_sq = anova.loc['C(domain)', 'sum_sq'] / anova['sum_sq'].sum()
print(f"η² = {eta_sq:.3f}")

def random_score_var(n_genes, n_iter=1000):
    vars_ = []
    for _ in range(n_iter):
        rand_genes = random.sample(list(adata.var_names), n_genes)
        X = adata[:, rand_genes].X
        if sparse.issparse(X): X = X.toarray()
        med = np.median(X, 0); mad = np.median(np.abs(X - med), 0); mad[mad==0]=1e-12
        vars_.append(((X - med) / mad).mean(1).var())
    return np.array(vars_)

null_var = random_score_var(hits)
real_var = adata.obs[score_col].var()
p_val = (null_var >= real_var).mean()
print("p-value =", p_val)

marker = "MS4A1"
sc.pl.scatter(adata, x=score_col, y=marker)

lib_id = list(adata.uns["spatial"].keys())[0]
vmax99 = np.percentile(adata.obs[score_col], 99)
sq.pl.spatial_scatter(
    adata,
    library_id=lib_id,
    color=score_col,
    cmap="magma",
    img_alpha=0.4,
    size=1.2,
    vmax=vmax99,
    title=sig_name
)
plt.show()
