import pandas as pd

df = pd.read_parquet(
    r"F:\24-25\Final_project\Visium_HD_Human_Lung_Cancer_binned_outputs\binned_outputs\square_008um\spatial\tissue_positions.parquet"
)

columns_order = [
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres"
]
df = df[columns_order]

df.to_csv(
    r"F:\24-25\Final_project\Visium_HD_Human_Lung_Cancer_binned_outputs\binned_outputs\square_008um\spatial\tissue_positions_list.csv",
    index=False,
    header=False
)

