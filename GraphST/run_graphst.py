import scanpy as sc
import torch
import pandas as pd
from GraphST.GraphST import GraphST
from GraphST.utils import refine_label
import matplotlib.pyplot as plt
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')


from sklearn.decomposition import PCA
import numpy as np



tool = 'leiden'
radius = 50
n_clusters = 7

adata = sc.read_h5ad(r"F:\24-25\Final_project\OneDrive_2_2025-5-29\P055_SpatialTranscriptomics\Data\SpatialProject\GraphST-main\subset_14854X18085.h5ad")

from GraphST.preprocess import preprocess, construct_interaction_KNN

preprocess(adata)
adata = adata[:, adata.var['highly_variable']]
adata.var_names_make_unique()
adata.X = adata.X.astype('float32')
construct_interaction_KNN(adata, n_neighbors=10)


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = GraphST(
    adata=adata,
    device=device,
    learning_rate=0.001,
    epochs=600,
    dim_input=adata.shape[1],
    dim_output=64,
    datatype="10X",
)

adata_out = model.train()

from GraphST.utils import clustering

clustering(adata_out, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=True)
adata_out.write("subset_14854X18085_leiden_preprocessing_pca_7clu_dim64_radius80_ngbrs10_graphst.h5ad")


