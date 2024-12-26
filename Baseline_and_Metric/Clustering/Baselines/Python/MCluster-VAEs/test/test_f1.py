from time import perf_counter

import numpy as np
import pandas as pd
import seaborn as sns

from src.evaluation import cluster_f1, cluster_f1_


res = []
for nc in (list(range(2, 11)) + [15, 20, 30]):
    for ns in [i*100 for i in range(1, 10)]:
        print("nc=%d, ns=%d" % (nc, ns))
        for i in range(100):
            a = np.random.randint(0, nc, size=(ns,))
            b = np.random.randint(0, nc, size=(ns,))
            t1 = perf_counter()
            f1_1 = cluster_f1(a, b)
            t2 = perf_counter()
            tt1 = t2 - t1
            t1 = perf_counter()
            f1_2 = cluster_f1_(a, b)
            t2 = perf_counter()
            tt2 = t2 - t1
            assert f1_1 == f1_2
            # tqdm.write("%.4f, %.4f" % (f1_1, f1_2))
            res.append(dict(nc=nc, ns=ns, t1=tt1, t2=tt2, re=i))

df = pd.DataFrame.from_records(res)
df = df.melt(id_vars=["nc", "ns", "re"], var_name="func")
fg = sns.relplot(data=df, x="ns", y="value",
                 hue="func", kind="line", col="nc", col_wrap=3)
fg.savefig("./test_f1.png")
