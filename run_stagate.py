import scanpy as sc
import STAGATE
from STAGATE.Train_STAGATE import train_STAGATE
import numpy as np
import pandas as pd
import pickle
import tensorflow.compat.v1 as tf
tf.disable_eager_execution()


adata = sc.read_h5ad(r"F:\24-25\Final_project\OneDrive_2_2025-5-29\P055_SpatialTranscriptomics\Data\SpatialProject\GraphST-main\subset_14854X18085.h5ad")
adata.var_names_make_unique()


sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.raw = adata.copy()

adata = adata[:, adata.var['highly_variable']]

STAGATE.Cal_Spatial_Net(adata, rad_cutoff=30)

adata = STAGATE.train_STAGATE(
    adata,
    hidden_dims=[512, 32],
    alpha=0.6,
    n_epochs=1000,
    lr=0.001,
    key_added='STAGATE',
    save_attention=False,
    save_loss=True,
    save_reconstrction=True
)

sc.pp.neighbors(adata, use_rep="STAGATE")

for res in np.arange(0.4, 0.9, 0.05):
    sc.tl.leiden(adata, resolution=round(res, 3), key_added=f"leiden_{round(res, 3)}")
    print(f"res={round(res,3)} → n_clusters={adata.obs[f'leiden_{round(res,3)}'].nunique()}")


sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.65, key_added="domain")

adata.write("subset_14854X18085_7clu_alpha6_rad30_epoch1000_dim32_stagate.h5ad")





# loss_curve = adata.uns["STAGATE_loss"]
#
# plt.figure(figsize=(6,4))
# plt.plot(loss_curve, label="STAGATE training loss")
# plt.xlabel("Epoch")
# plt.ylabel("Loss")
# plt.title("STAGATE reconstruction loss curve")
# plt.legend()
# plt.show()
