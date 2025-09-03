import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
import matplotlib
matplotlib.use('TkAgg')

file_fold = r"F:\24-25\Final_project\Visium_HD_Human_Lung_Cancer_binned_outputs\binned_outputs\square_008um"

adata = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True)
x_min, x_max = 8500, 9600
y_min, y_max = 6600, 7999

subset = adata[
    (adata.obsm['spatial'][:, 0] >= x_min) & (adata.obsm['spatial'][:, 0] <= x_max) &
    (adata.obsm['spatial'][:, 1] >= y_min) & (adata.obsm['spatial'][:, 1] <= y_max)
].copy()

print(subset)

save_path = r"F:\24-25\Final_project\OneDrive_2_2025-5-29\P055_SpatialTranscriptomics\Data\SpatialProject\GraphST-main\subset_14854X18085_ROI1.h5ad"
subset.write(save_path)

# spatial = adata.obsm['spatial']
#
# fig, ax = plt.subplots(figsize=(10, 10))
# ax.scatter(spatial[:, 0], spatial[:, 1], s=2, color='steelblue')
# ax.set_title("Full spatial layout - decide ROI", fontsize=14)
# ax.set_xlabel("X")
# ax.set_ylabel("Y")
# ax.invert_yaxis()
# ax.grid(True)
#
# x_min, x_max = 7000, 10000
# y_min, y_max = 5000, 8000
# roi_rect = patches.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min,
#                              linewidth=2, edgecolor='red', facecolor='none')
# ax.add_patch(roi_rect)
#
# plt.tight_layout()
# plt.show()
# library_id = list(adata.uns["spatial"].keys())[0]
# hires_img = adata.uns["spatial"][library_id]["images"]["hires"]
# scale = adata.uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"]
#
# x_min, x_max = 7000, 10000
# y_min, y_max = 5000, 8000
#
# x0 = int(x_min * scale)
# x1 = int(x_max * scale)
# y0 = int(y_min * scale)
# y1 = int(y_max * scale)
#
# roi_img = hires_img[y0:y1, x0:x1, :]
# plt.figure(figsize=(6, 6))
# plt.imshow(roi_img)
# plt.axis("off")
# plt.tight_layout(pad=0)
# plt.show()

