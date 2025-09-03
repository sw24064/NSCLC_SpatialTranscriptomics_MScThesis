import scanpy as sc
import pandas as pd

adata_raw = sc.read_h5ad(r"F:\24-25\Final_project\OneDrive_2_2025-5-29\P055_SpatialTranscriptomics\Data\SpatialProject\GraphST-main\subset_14854X18085.h5ad")

adata_trained = sc.read_h5ad("subset_14854X18085_6clu_alpha6_stagate.h5ad")

adata_raw = adata_raw[adata_trained.obs_names, :]

# .T
expr_df = pd.DataFrame(
    adata_raw.X.T.toarray(),
    index=adata_raw.var_names,
    columns=adata_raw.obs_names
)
expr_df.to_csv("expression_matrix_stagate.csv")

