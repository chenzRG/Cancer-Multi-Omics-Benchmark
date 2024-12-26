# Multi-omics Cluster Variational Autoencoders (MCluster-VAEs)

In this study, we propose Multi-omics Clustering Variational Autoencoders (**MCluster-VAEs**), a novel end-to-end deep learning-based method, to perform clustering using **multi-omics** data. MCluster-VAEs describes multi-omics data using a new probabilistic model with a **global discrete latent variable** considered as the clusters. Using **Variational Bayes** approach, we derive a unified single-model architecture and a novel objective function, to fit the multi-omics data and infer the posterior probability of the clustering assignments. The results on TCGA database with four kinds of omics demonstrate that MCluster-VAEs performs significantly better than 10 other methods. In particular, MCluster-VAEs achieves an **unsupervised clustering accuracy** of about **0.82** and **adjusted Rand index** of about **0.78** on the data with known labels, an increase of 6%-75% compared to other methods.

![](./imgs/figure1.png)

Overview of MCluster-VAEs. (A) Two deep learning-based multi-omics clustering pipelines. (B) The data generative model. (C) The posterior variational inference path. (D) The architecture of MCluster-VAEs model for the experiments using mRNA, miRNA, CNA and methylation.

## Dependencies

1. Download and install [miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html).

2. Create an python environment and install required packages.

```bash
conda create -n mcluster && conda activate mcluster
bash environment.sh
```
> The environment can also be created by `create env -f enviroment.yml`

## Running

1. Prepare data. Example data (Pan Cancer) is in *release* and can be download to test.

2. Modify the configuration. The all configurations are in `conf/`, are stored in *yaml* format and managed by `hydra` module. The modifying part is the data path and output path, like following:

    ```yaml
    # config.yaml
    hydra:
        ...
        run:
            # indicate the root path contained outputs
            dir: output_path/
    dataset:
        # indicate the root path contained datasets
        root: data_path/
        ...
    ```

    The more configuration can be found in `conf/`, including data pipeline (`conf/dataset/`), model architecture (`conf/model/`) and training setttings (`conf/train/`). These configuration can be added in `conf/config.yaml`, too. The rules of *hydra configuration* can be found in [hydra documents](https://hydra.cc/).


3. Run script.

    ```bash
    python main.py
    ```

    The results will be in `output_path/`


## contact

rongzhiwei@bjmu.edu.cn or 18845728185@163.com
