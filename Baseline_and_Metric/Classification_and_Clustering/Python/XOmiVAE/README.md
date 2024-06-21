# XOmiVAE
[![GitHub license](https://img.shields.io/github/license/Naereen/StrapDown.js.svg)](https://github.com/zhangxiaoyu11/XOmiVAE/blob/main/LICENSE)
![Safe](https://img.shields.io/badge/Stay-Safe-red?logo=data:image/svg%2bxml;base64,PHN2ZyBpZD0iTGF5ZXJfMSIgZW5hYmxlLWJhY2tncm91bmQ9Im5ldyAwIDAgNTEwIDUxMCIgaGVpZ2h0PSI1MTIiIHZpZXdCb3g9IjAgMCA1MTAgNTEwIiB3aWR0aD0iNTEyIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciPjxnPjxnPjxwYXRoIGQ9Im0xNzQuNjEgMzAwYy0yMC41OCAwLTQwLjU2IDYuOTUtNTYuNjkgMTkuNzJsLTExMC4wOSA4NS43OTd2MTA0LjQ4M2g1My41MjlsNzYuNDcxLTY1aDEyNi44MnYtMTQ1eiIgZmlsbD0iI2ZmZGRjZSIvPjwvZz48cGF0aCBkPSJtNTAyLjE3IDI4NC43MmMwIDguOTUtMy42IDE3Ljg5LTEwLjc4IDI0LjQ2bC0xNDguNTYgMTM1LjgyaC03OC4xOHYtODVoNjguMThsMTE0LjM0LTEwMC4yMWMxMi44Mi0xMS4yMyAzMi4wNi0xMC45MiA0NC41LjczIDcgNi41NSAxMC41IDE1LjM4IDEwLjUgMjQuMnoiIGZpbGw9IiNmZmNjYmQiLz48cGF0aCBkPSJtMzMyLjgzIDM0OS42M3YxMC4zN2gtNjguMTh2LTYwaDE4LjU1YzI3LjQxIDAgNDkuNjMgMjIuMjIgNDkuNjMgNDkuNjN6IiBmaWxsPSIjZmZjY2JkIi8+PHBhdGggZD0ibTM5OS44IDc3LjN2OC4wMWMwIDIwLjY1LTguMDQgNDAuMDctMjIuNjQgNTQuNjdsLTExMi41MSAxMTIuNTF2LTIyNi42NmwzLjE4LTMuMTljMTQuNi0xNC42IDM0LjAyLTIyLjY0IDU0LjY3LTIyLjY0IDQyLjYyIDAgNzcuMyAzNC42OCA3Ny4zIDc3LjN6IiBmaWxsPSIjZDAwMDUwIi8+PHBhdGggZD0ibTI2NC42NSAyNS44M3YyMjYuNjZsLTExMi41MS0xMTIuNTFjLTE0LjYtMTQuNi0yMi42NC0zNC4wMi0yMi42NC01NC42N3YtOC4wMWMwLTQyLjYyIDM0LjY4LTc3LjMgNzcuMy03Ny4zIDIwLjY1IDAgNDAuMDYgOC4wNCA1NC42NiAyMi42NHoiIGZpbGw9IiNmZjRhNGEiLz48cGF0aCBkPSJtMjEyLjgzIDM2MC4xMnYzMGg1MS44MnYtMzB6IiBmaWxsPSIjZmZjY2JkIi8+PHBhdGggZD0ibTI2NC42NSAzNjAuMTJ2MzBoMzYuMTRsMzIuMDQtMzB6IiBmaWxsPSIjZmZiZGE5Ii8+PC9nPjwvc3ZnPg==)
[![GitHub stars](https://img.shields.io/github/stars/zhangxiaoyu11/XOmiVAE.svg?style=social&label=Star&maxAge=2592000)](https://github.com/zhangxiaoyu11/XOmiVAE/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/zhangxiaoyu11/XOmiVAE.svg?style=social&label=Fork&maxAge=2592000)](https://github.com/zhangxiaoyu11/XOmiVAE/network/members)

**XOmiVAE: An Interpretable Deep Learning Framework for Cancer Classification using High-dimensional Omics Data**

**Xiaoyu Zhang** (x.zhang18@imperial.ac.uk)

Data Science Institute, Imperial College London

### Introduction
-  XOmiVAE is a novel interpretable deep learning model for cancer classification using high-dimensional omics data.
-  XOmiVAE provides contribution of each input molecular feature and latent dimension to the prediction.
-  XOmiVAE is able to explain unsupervised clusters produced by the VAE clustering part of the network.
-  XOmiVAE explanations of the downstream prediction were evaluated by biological annotation and literature, which aligned with current domain knowledge.
-  XOmiVAE shows great potential for novel biomedical knowledge discovery from deep learning models.

For examples of how to implement the code, see main.py. The explainable functions are in omiShapExplainer.py.

Paper Link: [Briefings in Bioinformatics](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbab315/6353242)

## Citation
If you use this code in your research, please cite our paper.
```bibtex
@Article{XOmiVAE2021,
    AUTHOR = {Withnell, Eloise and Zhang, Xiaoyu and Sun, Kai and Guo, Yike},
    TITLE = {XOmiVAE: an interpretable deep learning model for cancer classification using high-dimensional omics data},
    JOURNAL = {Briefings in Bioinformatics},
    YEAR = {2021},
    ISSN = {1477-4054},
    DOI = {10.1093/bib/bbab315},
    URL = {https://doi.org/10.1093/bib/bbab315},
}
```

## OmiEmbed
***Please also check our multi-task deep learning framework for multi-omics data***: 
[OmiEmbed](https://github.com/zhangxiaoyu11/OmiEmbed)

## License
This source code is licensed under the [MIT](https://github.com/zhangxiaoyu11/XOmiVAE/blob/main/LICENSE) license.
