The experimental data come from https://xenabrowser.net/datapages/?cohort=Visium%20HD%20Human%20Lung%20Cancer&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

Step 1
- Select a region using subset.py

Step 2
- STAGATE training was performed using https://github.com/zhanglabtools/STAGATE in combination with run_stagate.py
- GraphST training was performed using https://github.com/JinmiaoChenLab/GraphST in combination with run_graphst.py
- Use visual_follow.py to perform a preliminary visualization of the trained data. (Reference: GraphST 10X Visium — https://deepst‑tutorials.readthedocs.io/en/latest/Tutorial%201_10X%20Visium.html）(Reference: STAGATE 10X Visium — https://stagate.readthedocs.io/en/latest/T1_DLPFC.html)

Step 3
- Use check_result.py to evaluate the trained data.
- Use markergene.py to extract significant genes.
- Use gene_visual.py to visualize the extracted results.
- Use ROI.py to perform localized analysis on the regions of interest.
- Use GO.py to perform GO and GSEA enrichment analysis.
- Use further_analyz.py to perform in-depth analysis on specific mixed clusters, such as reclustering and adjusting the resolution.
- Use singleR.py to visualize the results of SingleR (SingleR is implemented in R and requires R version 4.3).
- Use signature.py to perform signature analysis (Using signature.csv).
- Use check_signature.py to analyze and visualize the signature results.

The PNG folder contains experimental result figures, which serve as supplementary results that were not included in the paper due to page limitations.

Other data in Onedrive:https://uob-my.sharepoint.com/:u:/g/personal/sw24064_bristol_ac_uk/ETyERgfh84xDrkjKj5Bg4skBCrkGRVZKLbahiv6Y_4pVqg?e=mslveg
Include:
The cropped h5ad file selected for the experiment.
The marker gene, GO, and GSEA enrichment result files were used for visualization.
The input file for SingleR.
One GraphST-trained h5ad file and one STAGATE-trained h5ad file were used for downstream analysis.
