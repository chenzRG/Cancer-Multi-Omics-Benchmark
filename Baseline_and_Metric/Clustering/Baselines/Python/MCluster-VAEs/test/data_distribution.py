import os

import pandas as pd
import seaborn as sns


data_dir = "/mnt/data1/share_data/TCGA/data_for_cluster/pancancer/fea/"
# clin_dir = "/mnt/data1/share_data/TCGA/data_for_cluster/pancancer/clinic/"

dfall = []
for n in ["meth", "rna", "CN", "miRNA"]:
    fn = os.path.join(data_dir, "PAN", n+".fea")
    df = pd.read_csv(fn, index_col=0).T
    df = df.sample(n=5, axis=1)
    df.columns = range(5)
    df["omic"] = n
    dfall.append(df)
dfall = pd.concat(dfall, axis=0)
dfall = dfall.melt(id_vars=["omic"], var_name="individual")

fg = sns.displot(dfall, x="value", hue="individual",
                 col="omic", col_wrap=2, kind="kde")
fg.savefig("./test/data_dist.png")
