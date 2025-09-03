import scanpy as sc
import numpy as np
from scipy import sparse
from scipy.sparse import csgraph
from sklearn.neighbors import NearestNeighbors
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 240)
pd.set_option('display.expand_frame_repr', False)


FILES_GRAPHST = {
    "neighbours10": "F:\\24-25\\Final_project\\OneDrive_2_2025-5-29\\P055_SpatialTranscriptomics\\Data\\GraphST-main\\GraphST-main\\subset_14854X18085_leiden_preprocessing_pca_7clu_graphst.h5ad",
    "neighbours12": "F:\\24-25\\Final_project\\OneDrive_2_2025-5-29\\P055_SpatialTranscriptomics\\Data\\GraphST-main\\GraphST-main\\subset_14854X18085_leiden_preprocessing_pca_7clu_dim64_radius50_ngbrs12_graphst.h5ad",
    "neighbours15": "F:\\24-25\\Final_project\\OneDrive_2_2025-5-29\\P055_SpatialTranscriptomics\\Data\\GraphST-main\\GraphST-main\\subset_14854X18085_leiden_preprocessing_pca_7clu_dim64_radius50_ngbrs15_graphst.h5ad",
}
# FILES_GRAPHST = {
#     "cluster7": "F:\\24-25\\Final_project\\OneDrive_2_2025-5-29\\P055_SpatialTranscriptomics\\Data\\GraphST-main\\GraphST-main\\subset_14854X18085_leiden_preprocessing_pca_7clu_graphst.h5ad",
#     "cluster8": "F:\\24-25\\Final_project\\OneDrive_2_2025-5-29\\P055_SpatialTranscriptomics\\Data\\GraphST-main\\GraphST-main\\subset_14854X18085_leiden_preprocessing_pca_8clu_graphst.h5ad",
#     "cluster9": "F:\\24-25\\Final_project\\OneDrive_2_2025-5-29\\P055_SpatialTranscriptomics\\Data\\GraphST-main\\GraphST-main\\subset_14854X18085_leiden_preprocessing_pca_9clu_graphst.h5ad",
# }

CLUSTER_KEY = "leiden"
EMBED_KEY   = "emb_pca"
K_LIST      = [6, 7, 8]


def _spatial_knn_indices(coords: np.ndarray, k: int):
    nbrs = NearestNeighbors(n_neighbors=k + 1, metric="euclidean")
    nbrs.fit(coords)
    ind = nbrs.kneighbors(return_distance=False)[:, 1:]
    return ind

def neighbor_purity(adata, labels, k=8):
    coords = np.asarray(adata.obsm["spatial"])
    ind = _spatial_knn_indices(coords, k)
    lab = np.asarray(labels)
    same = (lab[ind] == lab[:, None])
    return float(same.mean())

def boundary_ratio(adata, labels, k=8):
    coords = np.asarray(adata.obsm["spatial"])
    ind = _spatial_knn_indices(coords, k)
    n = adata.n_obs
    rows = np.repeat(np.arange(n), k)
    cols = ind.ravel()
    m = rows < cols
    rows, cols = rows[m], cols[m]
    lab = np.asarray(labels)
    diff = (lab[rows] != lab[cols])
    return float(diff.sum() / diff.size)

def connectivity_stats(adata, labels, k=8):
    coords = np.asarray(adata.obsm["spatial"])
    ind = _spatial_knn_indices(coords, k)
    n = adata.n_obs
    rows = np.repeat(np.arange(n), k)
    cols = ind.ravel()
    data = np.ones_like(rows, dtype=np.int8)
    A = sparse.coo_matrix((data, (rows, cols)), shape=(n, n))
    A = A.maximum(A.T).tocsr()

    labs = pd.Categorical(labels)
    lab_arr = labs.codes
    cats = labs.categories

    per_cluster = {}
    weights = []
    values = []
    for i, cat in enumerate(cats):
        mask = (lab_arr == i)
        size = int(mask.sum())
        if size <= 1:
            per_cluster[str(cat)] = 1.0
            weights.append(size)
            values.append(1.0)
            continue
        sub = A[mask][:, mask]
        n_comp, comp = csgraph.connected_components(sub, directed=False)
        sizes = np.bincount(comp)
        if sizes.size == 1:
            ratio = 1.0
        else:
            top2 = np.sort(sizes)[-2:].sum()
            ratio = float(top2 / sizes.sum())
        per_cluster[str(cat)] = ratio
        weights.append(size)
        values.append(ratio)

    weighted_mean = float(np.average(values, weights=np.clip(weights, 1, None)))
    unweighted_mean = float(np.mean(values))
    pct_ge_0_90 = float(np.mean([v >= 0.90 for v in values]))
    pct_ge_0_80 = float(np.mean([v >= 0.80 for v in values]))
    return per_cluster, weighted_mean, unweighted_mean, pct_ge_0_90, pct_ge_0_80

def silhouette_in_embedding(adata, labels, embed_key="STAGATE", n_comps=50):
    if embed_key in adata.obsm:
        X = adata.obsm[embed_key]
    elif "X_pca" in adata.obsm:
        X = adata.obsm["X_pca"]
    else:
        X = sc.pp.pca(adata, n_comps=n_comps, copy=True).obsm["X_pca"]
    y = pd.Categorical(labels).codes
    if len(np.unique(y)) < 2:
        return np.nan
    return float(silhouette_score(X, y, metric="euclidean"))

def ari_nmi_between(adata_a, key_a, adata_b, key_b):
    common = adata_a.obs_names.intersection(adata_b.obs_names)
    if len(common) == 0:
        raise ValueError("cannot compare")
    la = adata_a[common].obs[key_a].astype(str).values
    lb = adata_b[common].obs[key_b].astype(str).values
    ari = float(adjusted_rand_score(la, lb))
    nmi = float(normalized_mutual_info_score(la, lb))
    return ari, nmi

def spot_stability_across_models(loaded_dict, cluster_key):
    tags = list(loaded_dict.keys())
    if len(tags) < 2:
        raise ValueError("2 models at least")
    common = None
    for ad in loaded_dict.values():
        common = ad.obs_names if common is None else common.intersection(ad.obs_names)
    common = common.sort_values()
    m = len(tags)
    n = len(common)
    labels = [loaded_dict[t][common].obs[cluster_key].astype(str).to_numpy() for t in tags]
    total_pairs = m * (m - 1) // 2
    eq_sum = np.zeros(n, dtype=float)
    for i in range(m):
        li = labels[i]
        for j in range(i + 1, m):
            eq_sum += (li == labels[j])
    stability = eq_sum / total_pairs
    s = pd.Series(stability, index=common, name="spot_stability")
    return s, 1.0 - s

def permutation_p_value(adata, labels, metric_fn, k=8, n_perm=200, random_state=0):
    rng = np.random.default_rng(random_state)
    stat = metric_fn(adata, labels, k=k)
    labs = np.asarray(labels).copy()
    null = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        rng.shuffle(labs)
        null[i] = metric_fn(adata, labs, k=k)
    mu = null.mean()
    sd = null.std(ddof=1) if null.std(ddof=1) > 1e-12 else 1e-12
    z = (stat - mu) / sd
    p = (np.sum(np.abs(null - mu) >= np.abs(stat - mu)) + 1.0) / (n_perm + 1.0)
    return float(stat), float(z), float(p)

def boundary_contrast(adata, labels, embed_key="emb_pca", k=8):
    if embed_key in adata.obsm:
        X = np.asarray(adata.obsm[embed_key])
    elif "X_pca" in adata.obsm:
        X = np.asarray(adata.obsm["X_pca"])
    else:
        X = sc.pp.pca(adata, n_comps=50, copy=True).obsm["X_pca"]

    coords = np.asarray(adata.obsm["spatial"])
    ind = _spatial_knn_indices(coords, k)
    n = adata.n_obs
    rows = np.repeat(np.arange(n), k)
    cols = ind.ravel()
    lab = np.asarray(labels)
    same = (lab[rows] == lab[cols])

    d = np.linalg.norm(X[rows] - X[cols], axis=1)
    d_same = d[same]
    d_cross = d[~same]
    return {
        "dist_same_mean": float(d_same.mean()) if d_same.size else np.nan,
        "dist_cross_mean": float(d_cross.mean()) if d_cross.size else np.nan,
        "contrast_ratio": float(d_cross.mean() / d_same.mean()) if d_same.size and d_cross.size else np.nan,
        "contrast_diff":  float(d_cross.mean() - d_same.mean()) if d_same.size and d_cross.size else np.nan,
    }


def _align_labels_to_ref(ref_labels: np.ndarray, other_labels: np.ndarray) -> np.ndarray:
    ref_cat = pd.Categorical(ref_labels)
    oth_cat = pd.Categorical(other_labels)
    R, O = ref_cat.categories.size, oth_cat.categories.size

    cm = np.zeros((R, O), dtype=int)
    r_codes, o_codes = ref_cat.codes, oth_cat.codes
    for i in range(ref_labels.size):
        cm[r_codes[i], o_codes[i]] += 1

    r_ind, c_ind = linear_sum_assignment(-cm)
    mapping = {c: r for r, c in zip(r_ind, c_ind)}

    mapped_codes = np.array([mapping.get(code, code) for code in o_codes], dtype=int)
    mapped_names = np.asarray(ref_cat.categories)[mapped_codes]
    return mapped_names.astype(str)

def spot_consensus_across_models_aligned(loaded_dict, cluster_key, ref_tag=None):
    tags = list(loaded_dict.keys())
    if len(tags) < 2:
        raise ValueError("2 models at least")
    if ref_tag is None:
        ref_tag = tags[0]

    common = None
    for ad in loaded_dict.values():
        common = ad.obs_names if common is None else common.intersection(ad.obs_names)
    common = common.sort_values()

    ref = loaded_dict[ref_tag][common].obs[cluster_key].astype(str).to_numpy()

    aligned = []
    for t in tags:
        lab = loaded_dict[t][common].obs[cluster_key].astype(str).to_numpy()
        if t == ref_tag:
            aligned.append(ref)
        else:
            aligned.append(_align_labels_to_ref(ref, lab))

    L = np.vstack(aligned)
    M = L.shape[0]
    consensus = []
    for j in range(L.shape[1]):
        _, cnt = np.unique(L[:, j], return_counts=True)
        consensus.append(cnt.max() / M)
    consensus = pd.Series(consensus, index=common, name="spot_consensus_aligned")
    return consensus, 1.0 - consensus


# main
rows = []
detail_connectivity = {}

loaded = {}

for tag, path in FILES_GRAPHST.items():
    ad = sc.read_h5ad(path)
    if "spatial" not in ad.obsm_keys():
        raise ValueError(f"{path} no obsm['spatial']")
    if CLUSTER_KEY not in ad.obs.columns:
        raise ValueError(f"{path} no obs['{CLUSTER_KEY}']")
    loaded[tag] = ad

for tag, adata in loaded.items():
    labels = adata.obs[CLUSTER_KEY].astype(str)

    purities, boundaries, conn_w_list, conn_u_list, pct90_list, pct80_list = [], [], [], [], [], []
    per_clu_any = None
    for k in K_LIST:
        purities.append(neighbor_purity(adata, labels, k=k))
        boundaries.append(boundary_ratio(adata, labels, k=k))
        per_clu, conn_w, conn_u, pct90, pct80 = connectivity_stats(adata, labels, k=k)
        per_clu_any = per_clu
        conn_w_list.append(conn_w); conn_u_list.append(conn_u)
        pct90_list.append(pct90);   pct80_list.append(pct80)

    detail_connectivity[tag] = per_clu_any

    sil = silhouette_in_embedding(adata, labels, embed_key=EMBED_KEY)

    rows.append({
        "model": tag,
        f"neighbor_purity@k{K_LIST}": float(np.mean(purities)),
        f"boundary_ratio@k{K_LIST}":  float(np.mean(boundaries)),
        "connectivity_weighted":      float(np.mean(conn_w_list)),
        "connectivity_unweighted":    float(np.mean(conn_u_list)),
        "pct_clusters_ratio>=0.90":   float(np.mean(pct90_list)),
        "pct_clusters_ratio>=0.80":   float(np.mean(pct80_list)),
        f"silhouette_on_{EMBED_KEY}": sil,
        "n_clusters": int(pd.Categorical(labels).categories.size),
        "n_obs": int(adata.n_obs),
    })

df = pd.DataFrame(rows).set_index("model").sort_index()
print("\n=== GraphST (k in {}) summary ===".format(K_LIST))
print(df.round(4).to_string())

print("\n=== Per-cluster connectivity (top-2 components / cluster size) ===")
for tag, d in detail_connectivity.items():
    print(f"\n[{tag}]")
    s = pd.Series(d, dtype=float).sort_index()
    print(s.round(4).to_string())

pairs = []
tags = list(loaded.keys())
for i in range(len(tags)):
    for j in range(i+1, len(tags)):
        A, B = loaded[tags[i]], loaded[tags[j]]
        common = A.obs_names.intersection(B.obs_names)
        la = A[common].obs[CLUSTER_KEY].astype(str).values
        lb = B[common].obs[CLUSTER_KEY].astype(str).values
        pairs.append({
            "A": tags[i], "B": tags[j],
            "ARI": adjusted_rand_score(la, lb),
            "NMI": normalized_mutual_info_score(la, lb),
        })
if pairs:
    print("\n=== Pairwise ARI/NMI across GraphST clusterings ===")
    print(pd.DataFrame(pairs).sort_values("ARI", ascending=False).round(4).to_string())

if len(loaded) >= 2:
    consensus_s, uncertainty_s = spot_consensus_across_models_aligned(loaded, CLUSTER_KEY, ref_tag=list(loaded.keys())[0])
    print("\n=== Spot-level consensus after label alignment (0-1) ===")
    print(consensus_s.describe().round(4).to_string())
    print("\n=== Spot-level uncertainty (1 - consensus) (0-1) ===")
    print(uncertainty_s.describe().round(4).to_string())

    k_test = K_LIST[-1] if len(K_LIST) > 0 else 8
    extra_rows = []
    for tag, ad in loaded.items():
        labs = ad.obs[CLUSTER_KEY].astype(str).values
        pur_stat, pur_z, pur_p = permutation_p_value(ad, labs, neighbor_purity, k=k_test, n_perm=200, random_state=0)
        br_stat, br_z, br_p = permutation_p_value(ad, labs, boundary_ratio, k=k_test, n_perm=200, random_state=1)
        bc = boundary_contrast(ad, labs, embed_key=EMBED_KEY, k=k_test)
        extra_rows.append({
            "model": tag,
            "k_for_tests": k_test,
            "purity_stat": pur_stat, "purity_z": pur_z, "purity_p": pur_p,
            "boundary_stat": br_stat, "boundary_z": br_z, "boundary_p": br_p,
            "dist_same_mean": bc["dist_same_mean"],
            "dist_cross_mean": bc["dist_cross_mean"],
            "boundary_contrast_ratio": bc["contrast_ratio"],
            "boundary_contrast_diff": bc["contrast_diff"],
        })

    extra_df = pd.DataFrame(extra_rows).set_index("model").sort_index()
    print("\n=== Spatial permutation tests & boundary contrast (k={}) ===".format(k_test))
    print(extra_df.round(4).to_string())


from sklearn.mixture import GaussianMixture

# ROI
def subset_by_roi(adata, x_min, x_max, y_min, y_max):
    S = np.asarray(adata.obsm["spatial"])
    keep = (S[:,0] >= x_min) & (S[:,0] <= x_max) & (S[:,1] >= y_min) & (S[:,1] <= y_max)
    return adata[keep].copy()

ROI = (8100, 9900, 6100, 7999)
FINE_TAG   = "neighbours10"
COARSE_TAG = "neighbours15"

loaded_roi = {}
for tag, ad in loaded.items():
    loaded_roi[tag] = subset_by_roi(ad, *ROI)
    print(f"[ROI] {tag}: n_spots = {loaded_roi[tag].n_obs}")

tags = list(loaded_roi.keys())
pairs = []
for i in range(len(tags)):
    for j in range(i+1, len(tags)):
        A, B = loaded_roi[tags[i]], loaded_roi[tags[j]]
        common = A.obs_names.intersection(B.obs_names)
        la = A[common].obs[CLUSTER_KEY].astype(str).values
        lb = B[common].obs[CLUSTER_KEY].astype(str).values
        pairs.append({"A": tags[i], "B": tags[j],
                      "ARI": adjusted_rand_score(la, lb),
                      "NMI": normalized_mutual_info_score(la, lb)})
if pairs:
    print("\n[ROI] Pairwise ARI/NMI across models")
    print(pd.DataFrame(pairs).sort_values("ARI", ascending=False).round(4).to_string(index=False))

def _get_latent(ad, key=EMBED_KEY):
    if key in ad.obsm:
        return np.asarray(ad.obsm[key])
    elif "X_pca" in ad.obsm:
        return np.asarray(ad.obsm["X_pca"])
    else:
        return sc.pp.pca(ad, n_comps=50, copy=True).obsm["X_pca"]

print("\n[ROI] Unsupervised evidence (GMM)")
for tag, ad in loaded_roi.items():
    X = _get_latent(ad, EMBED_KEY)
    gm1 = GaussianMixture(n_components=1, covariance_type="full", random_state=0).fit(X)
    gm2 = GaussianMixture(n_components=2, covariance_type="full", random_state=0).fit(X)
    bic1, bic2 = gm1.bic(X), gm2.bic(X)
    delta_bic = bic1 - bic2       # >10 强支持“两簇”
    # 用 GMM-2 的分配算 silhouette（仅参考）
    lab2 = gm2.predict(X)
    sil2 = silhouette_score(X, lab2) if len(np.unique(lab2)) > 1 else np.nan
    print(f"[{tag}] BIC(1)={bic1:.1f}, BIC(2)={bic2:.1f}, ΔBIC={delta_bic:.1f}, silhouette(k=2)={sil2:.3f}")

k_test = 8
rows = []
for tag, ad in loaded_roi.items():
    labs = ad.obs[CLUSTER_KEY].astype(str).values
    pur_stat, pur_z, pur_p = permutation_p_value(ad, labs, neighbor_purity, k=k_test, n_perm=200, random_state=0)
    br_stat,  br_z,  br_p  = permutation_p_value(ad, labs, boundary_ratio , k=k_test, n_perm=200, random_state=1)
    bc = boundary_contrast(ad, labs, embed_key=EMBED_KEY, k=k_test)
    rows.append({
        "model": tag, "k": k_test,
        "neighbor_purity": pur_stat, "purity_p": pur_p,
        "boundary_ratio": br_stat,  "boundary_p": br_p,
        "dist_same_mean": bc["dist_same_mean"],
        "dist_cross_mean": bc["dist_cross_mean"],
        "contrast_ratio": bc["contrast_ratio"],
        "contrast_diff":  bc["contrast_diff"],
        "n_clusters_in_ROI": int(pd.Categorical(labs).categories.size)
    })
print("\n[ROI] Spatial tests")
print(pd.DataFrame(rows).set_index("model").round(4).to_string())

def _topn_marker_set(ad, groupby, clabel, topn=50):
    sc.tl.rank_genes_groups(ad, groupby=groupby, method="wilcoxon")
    df = sc.get.rank_genes_groups_df(ad, group=clabel).sort_values("logfoldchanges", ascending=False)
    return set(df["names"].head(topn).values)

def _two_cluster_DE_counts(ad, groupby, cA, cB, logfc_thr=0.25, q_thr=0.05):
    sub = ad[ad.obs[groupby].isin([cA, cB])].copy()
    # 统一成 A vs B：把 B 当 reference
    sc.tl.rank_genes_groups(sub, groupby=groupby, groups=[cA], reference=cB, method="wilcoxon")
    df = sc.get.rank_genes_groups_df(sub, group=cA)
    sig = df[(df["pvals_adj"] < q_thr) & (np.abs(df["logfoldchanges"]) > logfc_thr)]
    return int(sig.shape[0])

fine_ad = loaded_roi[FINE_TAG]
counts = fine_ad.obs[CLUSTER_KEY].value_counts()
if counts.shape[0] >= 2:
    cA, cB = "2","6"
    n_sig = _two_cluster_DE_counts(fine_ad, CLUSTER_KEY, cA, cB, logfc_thr=0.25, q_thr=0.05)
    A_set = _topn_marker_set(fine_ad, CLUSTER_KEY, cA, topn=50)
    B_set = _topn_marker_set(fine_ad, CLUSTER_KEY, cB, topn=50)
    jacc_AB = len(A_set & B_set) / max(1, len(A_set | B_set))
    print(f"\n[ROI][{FINE_TAG}] Largest two clusters: {cA} vs {cB}")
    print(f" DE(sig) count (|logFC|>0.25, q<0.05): {n_sig}")
    print(f" Marker Jaccard (top50): {jacc_AB:.3f}")
else:
    cA = cB = None
    print(f"\n[ROI][{FINE_TAG}] Less than two clusters within ROI; skip DE/Jaccard.")

coarse_ad = loaded_roi[COARSE_TAG]
coarse_main = coarse_ad.obs[CLUSTER_KEY].value_counts().index[0]
C_set = _topn_marker_set(coarse_ad, CLUSTER_KEY, coarse_main, topn=50)

if cA is not None and cB is not None:
    jacc_A_C = len(A_set & C_set) / max(1, len(A_set | C_set))
    jacc_B_C = len(B_set & C_set) / max(1, len(B_set | C_set))
    print(f"[ROI] Marker Jaccard to coarse main (top50): A~C={jacc_A_C:.3f}, B~C={jacc_B_C:.3f}")

