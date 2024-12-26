<div align="center">
	<img src="hex-maui.png" alt="maui">
</div>

# maui

[![Downloads](https://pepy.tech/badge/maui-tools)](https://pepy.tech/project/maui-tools)  [![codecov](https://codecov.io/gh/bimsbbioinfo/maui/branch/master/graph/badge.svg)](https://codecov.io/gh/bimsbbioinfo/maui) [![Codacy Badge](https://api.codacy.com/project/badge/Grade/36c1f3f252b543139fd930ba5f674535)](https://www.codacy.com/app/jonathanronen/maui?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=BIMSBbioinfo/maui&amp;utm_campaign=Badge_Grade) [![PyPI version](https://badge.fury.io/py/maui-tools.svg)](https://badge.fury.io/py/maui-tools) [![Documentation Status](https://readthedocs.org/projects/maui/badge/?version=latest)](https://maui.readthedocs.io/en/latest/?badge=latest) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/python/black)



<!-- [![Build Status](https://travis-ci.com/BIMSBbioinfo/maui.svg?branch=master)](https://travis-ci.com/BIMSBbioinfo/maui) -->


Multi-omics Autoencoder Integration (**maui**) is a python package for multi-omics data analysis. It is based on a bayesian latent factor model, with inference done using artificial neural networks. For details, check out our LSA paper: https://www.life-science-alliance.org/content/2/6/e201900517

## Installation

maui works with Python 3.6 and TensorFlow 1.1 (does not yet support the yet unreleased TensorFlow 2.0). The easiest way to install is from pypi:

	pip install -U maui-tools

This will install all necessary dependencies including keras an tensorflow. The default tensorflow (cpu) will be installed. If tensorflow GPU is needed, please install it prior to installation of maui.

The development version may be installed by cloning this repo and running `python setup.py install`, or, using `pip` directly from github:

	pip install -e git+https://github.com/BIMSBbioinfo/maui.git#egg=maui


#### Optional dependencies

Survival analysis functionality supplied by lifelines [1]. It may be installed directly from pip using `pip install lifelines`.

## Usage

See [the vignette](vignette/maui_vignette.ipynb), and check out [the documentation](https://maui.readthedocs.io/en/latest/).


## Citation

>  Evaluation of colorectal cancer subtypes and cell lines using deep learning. Jonathan Ronen, Sikander Hayat, Altuna Akalin. Life Science Alliance Dec 2019, 2 (6) e201900517; DOI: 10.26508/lsa.201900517

## Contributing

Open an issue, send us a pull request, or shoot us an e-mail.

## License

maui is released under the [GNU General Public License v3.0](LICENSE) or later.

---------------------
@jonathanronen, BIMSBbioinfo, 2018


[1]: https://github.com/CamDavidsonPilon/lifelines
