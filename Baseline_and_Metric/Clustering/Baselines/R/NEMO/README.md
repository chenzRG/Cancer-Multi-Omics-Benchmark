# NEMO
This repository includes an implementation of the NEMO (NEighborhood based Multi-Omics clustering) multi-omics clustering algorithm.
The NEMO package is included in the NEMO directory. The other files, in the top directory, include code to reproduce all NEMO results in the paper. 
NEMO can be installed as follows:
```{r}
devtools::install_github('Shamir-Lab/NEMO/NEMO')
```
Or just download the files and run:
```{r}
source('NEMO.R')
```

Basic use of NEMO is as follows:

```{r}
# omics.list is a list of data frames, where in each data frame columns are samples and rows are features.
# note that colnames for each data frame should be set.
# Two toy datasets, omic1 and omic2, are also loaded whenever NEMO is loaded, but here we read them from a file.
omic1 = read.table(file='data/sample_omic1')
omic2 = read.table(file='data/sample_omic2')
omics.list = list(omic1, omic2)
clustering = nemo.clustering(omics.list)

# the number of clusters and number of nearest neighbors used can also be given as input.
# if they are not given, nemo decides the values itself.
nemo.clustering(omics.list, num.clusters=2, num.neighbors=50)
```
A more advanced use is as follows:
```{r}
# k can also be set to NA, in which case nemo chooses its value.
# nemo.affinity.graph is the integrated affinity graph.
affinity.graph = nemo.affinity.graph(omics.list, k = 20)

# ask nemo to estimate the number of clusters.
num.clusters = nemo.num.clusters(affinity.graph)

# clustering is the cluster assignment vector, ordered by the columns of affinity.graph.
clustering = spectralClustering(affinity.graph, num.clusters)
names(clustering) = colnames(affinity.graph)
```

NEMO requires prior installation of the R library SNFtool, and also uses some of its code:
https://cran.r-project.org/web/packages/SNFtool/index.html
