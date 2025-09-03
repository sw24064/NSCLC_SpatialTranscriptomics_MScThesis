import scanpy as sc
import numpy as np
from scipy import sparse
from scipy.sparse import csgraph
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score
import networkx as nx
from libpysal.weights import KNN
from esda.moran import Moran
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score
import pandas as pd
pd.set_option('display.max_columns', None)     # 列不折叠
pd.set_option('display.width', 240)            # 控制台宽一点
pd.set_option('display.expand_frame_repr', False)


FILES_STAGATE = {
    "alpha0.6": "F:\\24-25\\Final_project\\OneDrive_2_2025-5-29\\P055_SpatialTranscriptomics\\Data\\STAGATE-main\\STAGATE-main\\subset_14854X18085_7clu_alpha6_stagate.h5ad",
    "alpha0.7": "F:\\24-25\\Final_project\\OneDrive_2_2025-5-29\\P055_SpatialTranscriptomics\\Data\\STAGATE-main\\STAGATE-main\\subset_14854X18085_7clu_alpha7_stagate.h5ad",
    "alpha0.8": "F:\\24-25\\Final_project\\OneDrive_2_2025-5-29\\P055_SpatialTranscriptomics\\Data\\STAGATE-main\\STAGATE-main\\subset_14854X18085_7clu_alpha8_stagate.h5ad",
    "alpha0.9": "F:\\24-25\\Final_project\\OneDrive_2_2025-5-29\\P055_SpatialTranscriptomics\\Data\\STAGATE-main\\STAGATE-main\\subset_14854X18085_7clu_alpha9_stagate.h5ad",  # 你给的示例
    "alpha1.0": "F:\\24-25\\Final_project\\OneDrive_2_2025-5-29\\P055_SpatialTranscriptomics\\Data\\STAGATE-main\\STAGATE-main\\subset_14854X18085_7clu_alpha10_stagate.h5ad",
}
# FILES_STAGATE = {
#     "cluster6": "F:\\24-25\\Final_project\\OneDrive_2_2025-5-29\\P055_SpatialTranscriptomics\\Data\\STAGATE-main\\STAGATE-main\\subset_14854X18085_6clu_alpha6_stagate.h5ad",
#     "cluster7": "F:\\24-25\\Final_project\\OneDrive_2_2025-5-29\\P055_SpatialTranscriptomics\\Data\\STAGATE-main\\STAGATE-main\\subset_14854X18085_6clu_alpha8_stagate.h5ad",
#     "cluster9": "F:\\24-25\\Final_project\\OneDrive_2_2025-5-29\\P055_SpatialTranscriptomics\\Data\\STAGATE-main\\STAGATE-main\\subset_14854X18085_6clu_alpha9_stagate.h5ad",
# }

CLUSTER_KEY = "domain"
EMBED_KEY   = "STAGATE"
K_SPATIAL   = [6, 7, 8]


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

def modularity_score(adata, labels, k=8):
    coords = np.asarray(adata.obsm["spatial"])
    ind = _spatial_knn_indices(coords, k)
    n = adata.n_obs
    rows = np.repeat(np.arange(n), k)
    cols = ind.ravel()
    data = np.ones_like(rows, dtype=np.int8)
    A = sparse.coo_matrix((data, (rows, cols)), shape=(n, n))
    A = A.maximum(A.T).tocsr()
    G = nx.from_scipy_sparse_array(A)

    labs = pd.Categorical(labels)
    communities = [np.where(labs.codes == i)[0] for i in range(len(labs.categories))]
    return nx.algorithms.community.quality.modularity(G, communities)

def morans_I_score(adata, labels, k=8):
    coords = np.asarray(adata.obsm["spatial"])
    w = KNN.from_array(coords, k=k)
    y = pd.Categorical(labels).codes
    mi = Moran(y, w)
    return mi.I, mi.p_sim

# main
rows = []
detail_connectivity = {}

loaded = {}
for tag, path in FILES_STAGATE.items():
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
    for k in K_SPATIAL:
        purities.append(neighbor_purity(adata, labels, k=k))
        boundaries.append(boundary_ratio(adata, labels, k=k))
        per_clu, conn_w, conn_u, pct90, pct80 = connectivity_stats(adata, labels, k=k)
        per_clu_any = per_clu
        conn_w_list.append(conn_w); conn_u_list.append(conn_u)
        pct90_list.append(pct90);   pct80_list.append(pct80)

    detail_connectivity[tag] = per_clu_any

    sil = silhouette_in_embedding(adata, labels, embed_key=EMBED_KEY)
    X = adata.obsm[EMBED_KEY]
    y = pd.Categorical(labels).codes
    ch = calinski_harabasz_score(X, y) if len(np.unique(y)) > 1 else np.nan
    db = davies_bouldin_score(X, y) if len(np.unique(y)) > 1 else np.nan

    mod = np.mean([modularity_score(adata, labels, k=k) for k in K_SPATIAL])
    mi_val, mi_p = morans_I_score(adata, labels, k=np.max(K_SPATIAL))

    rows.append({
        "model": tag,
        f"neighbor_purity@k{K_SPATIAL}": float(np.mean(purities)),
        f"boundary_ratio@k{K_SPATIAL}": float(np.mean(boundaries)),
        "connectivity_weighted": float(np.mean(conn_w_list)),
        "connectivity_unweighted": float(np.mean(conn_u_list)),
        "pct_clusters_ratio>=0.90": float(np.mean(pct90_list)),
        "pct_clusters_ratio>=0.80": float(np.mean(pct80_list)),
        f"silhouette_on_{EMBED_KEY}": sil,
        "Calinski_Harabasz": ch,
        "Davies_Bouldin": db,
        "Modularity": mod,
        "Morans_I": mi_val,
        "Morans_p": mi_p,
        "n_clusters": int(pd.Categorical(labels).categories.size),
        "n_obs": int(adata.n_obs),
    })

df = pd.DataFrame(rows).set_index("model").sort_index()
print("\nSTAGATE (α sensitivity test) summary")
print(df.round(4).to_string())

print("\nPer-cluster connectivity")
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
    print("\nARI/NMI")
    print(pd.DataFrame(pairs).sort_values("ARI", ascending=False).round(4).to_string())