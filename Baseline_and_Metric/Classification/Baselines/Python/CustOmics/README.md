# CustOmics
[![GitHub license](https://img.shields.io/github/license/Naereen/StrapDown.js.svg)](https://github.com/HakimBenkirane/CustOMICS/LICENSE)
![Safe](https://img.shields.io/badge/Stay-Safe-red?logo=data:image/svg%2bxml;base64,PHN2ZyBpZD0iTGF5ZXJfMSIgZW5hYmxlLWJhY2tncm91bmQ9Im5ldyAwIDAgNTEwIDUxMCIgaGVpZ2h0PSI1MTIiIHZpZXdCb3g9IjAgMCA1MTAgNTEwIiB3aWR0aD0iNTEyIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciPjxnPjxnPjxwYXRoIGQ9Im0xNzQuNjEgMzAwYy0yMC41OCAwLTQwLjU2IDYuOTUtNTYuNjkgMTkuNzJsLTExMC4wOSA4NS43OTd2MTA0LjQ4M2g1My41MjlsNzYuNDcxLTY1aDEyNi44MnYtMTQ1eiIgZmlsbD0iI2ZmZGRjZSIvPjwvZz48cGF0aCBkPSJtNTAyLjE3IDI4NC43MmMwIDguOTUtMy42IDE3Ljg5LTEwLjc4IDI0LjQ2bC0xNDguNTYgMTM1LjgyaC03OC4xOHYtODVoNjguMThsMTE0LjM0LTEwMC4yMWMxMi44Mi0xMS4yMyAzMi4wNi0xMC45MiA0NC41LjczIDcgNi41NSAxMC41IDE1LjM4IDEwLjUgMjQuMnoiIGZpbGw9IiNmZmNjYmQiLz48cGF0aCBkPSJtMzMyLjgzIDM0OS42M3YxMC4zN2gtNjguMTh2LTYwaDE4LjU1YzI3LjQxIDAgNDkuNjMgMjIuMjIgNDkuNjMgNDkuNjN6IiBmaWxsPSIjZmZjY2JkIi8+PHBhdGggZD0ibTM5OS44IDc3LjN2OC4wMWMwIDIwLjY1LTguMDQgNDAuMDctMjIuNjQgNTQuNjdsLTExMi41MSAxMTIuNTF2LTIyNi42NmwzLjE4LTMuMTljMTQuNi0xNC42IDM0LjAyLTIyLjY0IDU0LjY3LTIyLjY0IDQyLjYyIDAgNzcuMyAzNC42OCA3Ny4zIDc3LjN6IiBmaWxsPSIjZDAwMDUwIi8+PHBhdGggZD0ibTI2NC42NSAyNS44M3YyMjYuNjZsLTExMi41MS0xMTIuNTFjLTE0LjYtMTQuNi0yMi42NC0zNC4wMi0yMi42NC01NC42N3YtOC4wMWMwLTQyLjYyIDM0LjY4LTc3LjMgNzcuMy03Ny4zIDIwLjY1IDAgNDAuMDYgOC4wNCA1NC42NiAyMi42NHoiIGZpbGw9IiNmZjRhNGEiLz48cGF0aCBkPSJtMjEyLjgzIDM2MC4xMnYzMGg1MS44MnYtMzB6IiBmaWxsPSIjZmZjY2JkIi8+PHBhdGggZD0ibTI2NC42NSAzNjAuMTJ2MzBoMzYuMTRsMzIuMDQtMzB6IiBmaWxsPSIjZmZiZGE5Ii8+PC9nPjwvc3ZnPg==)
[![GitHub stars](https://img.shields.io/github/stars/HakimBenkirane/CustOmics.svg?style=social&label=Star&maxAge=2592000)](https://github.com/HakimBenkirane/CustOmics/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/HakimBenkirane/CustOmics.svg?style=social&label=Fork&maxAge=2592000)](https://github.com/HakimBenkirane/CustOmics/network/members)

**CustOmics: A versatile deep-learning based strategy for multi-omics integration**

**Hakim Benkirane** (hakim.benkirane@centralesupelec.fr)

Oncostat Team, U1018 Inserm, CESP
Laboratory of mathematics and informatics of CentraleSupelec


### Introduction
-  CustOmics is a novel architecture for classification and survival outcome prediction.
-  CustOmics uses a new integration strategy for a more versatile multi-omics integration.
-  CustOmics is able to provide both end-to-end prediction and unsupervised latent representation.
-  CustOmics has been evaluated using multiple test cases for classification and survival using TCGA datasets.
-  CustOmics is able to explain, to a certain degree, classification results.


Paper Link: [Link to the published paper]()


### Downloading TCGA Data
To download omics data (formatted as .tsv files) and other clinical metadata, please refer to the [NIH Genomic Data Commons Data Portal](https://portal.gdc.cancer.gov/) and the [cBioPortal](https://www.cbioportal.org/).

### Running Experiments

Experiments can be executed through the script **main.py**, the basic usage to run a tumor type classification on the Pancancer dataset is as follows:

``` shell
python main.py --cohorts PANCAN --sources CNV,RNAseq,methyl --task classification --data_directory DATA_DIRECTORY --result_directory RESULTS_DIRECTORY
```

To run PAM50 classification task on TCGA-BRCA dataset:

``` shell
python main.py --cohorts TCGA-BRCA --sources CNV,RNAseq,methyl --task classification --data_directory DATA_DIRECTORY --result_directory RESULTS_DIRECTORY
```

To run survival tasks on specific datasets:

``` shell
python main.py --cohorts TCGA-BLCA,TCGA-BRCA,TCGA-LUAD,TCGA-GBM,TCGA-UCEC --sources CNV,RNAseq,methyl --task survival --data_directory DATA_DIRECTORY --result_directory RESULTS_DIRECTORY
```


## License
This source code is licensed under the [MIT](https://github.com/HakimBenkirane/CustOmics/LICENSE) license.

# Cite us


## Citation
If you use this code in your research, please cite our paper.
```bibtex
@article{benkirane2023,
    doi = {10.1371/journal.pcbi.1010921},
    author = {Benkirane, Hakim AND Pradat, Yoann AND Michiels, Stefan AND Cournède, Paul-Henry},
    journal = {PLOS Computational Biology},
    publisher = {Public Library of Science},
    title = {CustOmics: A versatile deep-learning based strategy for multi-omics integration},
    year = {2023},
    month = {03},
    volume = {19},
    url = {https://doi.org/10.1371/journal.pcbi.1010921},
    pages = {1-19},
    number = {3}
}
```
