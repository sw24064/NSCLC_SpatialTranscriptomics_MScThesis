# GO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('TkAgg')

file_path = "cluster0_subclustered_2clu/GO_cluster_0/GO_BP_enrichr.csv"
cluster_id = "0"

df = pd.read_csv(file_path)
df = df[df["Adjusted P-value"] < 0.05]
df = df[df["Combined Score"] > 30]
df = df.sort_values("Combined Score", ascending=False).head(10)
# df = df.reindex(df["NES"].abs().sort_values(ascending=False).index)
# df_top = df.head(10)

plt.figure(figsize=(8, 6))
sns.barplot(data=df, x="Combined Score", y="Term", palette="magma")

plt.title(f"Top 10 GO Terms for Subluster {cluster_id} in Cluster 0 (GraphST)")
plt.xlabel("Combined Score")
plt.ylabel("GO Term")
plt.tight_layout()
plt.show()
