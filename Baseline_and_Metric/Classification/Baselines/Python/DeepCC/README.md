# DeepCC
DeepCC: a novel deep learning-based framework for cancer molecular subtype classification

([See Gao, et al. Oncogenesis, 2019](https://www.nature.com/articles/s41389-019-0157-8))

### Updata
- Added supporting of deep learning with keras framework.
- Modified the deep netwrok architecture for better performance.

### Dependencies
The current version of DeepCC dependes on keras framework, with CPU and GPU supporting. Please install [keras](https://keras.rstudio.com/) in R first accroding to the instruction and make sure that keras able to run correctly.

### Installation
You can install DeepCC from GitHub directly using devtools.
```
install.packages("devtools")
devtools::install_github("CityUHK-CompBio/DeepCC")
```
or 
```
install.packages("remote")
remotes::install_github("CityUHK-CompBio/DeepCC")
```
### Quick start
As a case study, you can obtain well organized colorectal cancer data from CRCSC's repository on [Synapse](https://www.synapse.org/#!Synapse:syn2623706/wiki/).

DeepCC only need two input for start.
- a data.frame containing gene expression profiles (each row represents one patient sample and the column names should be Entrez identifiers of genes)
- a character vector containing training labels (`NA` is allowed since DeepCC can ignore them automatically)

Now assume you have the gene expression profiles in `eps` and training labels in `labels`.
```
library(DeepCC)
library(keras)

# get functional spectra from gene expression profiles and the geneSets will use MSigDBv7 by default. Add parameter `geneSets = "MSigDBvx"` (x = 5/6/7) or `geneSets = MSigDBr` (MSigDBr attained from R package `msigdbr` with function get_msigdbr() in DeepCC) for different version of MSigDB.
fs <- getFunctionalSpectra(eps)

# train DeepCC model
deepcc_model <- train_DeepCC_model(fs, labels)

# obtain deep features 
df <- get_DeepCC_features(deepcc_model, fs)
```

After training, now you can use your DeepCC model to classify new sample(s). DeepCC can classify samples in a data set, as well as individual samples. The input data should be in the same format as above gene expression profile(s).

```
# classify new data set use trained DeepCC model
# for a batch of samples
new_fs <- getFunctionalSpectra(new_eps)
pred_labels <- get_DeepCC_label(deepcc_model, new_fs)

# for a given single sample, you have to provide a reference expression profile.
new_fs_single <- getFunctionalSpectrum(new_ep, refExp = "COADREAD")
pred_label <- get_DeepCC_label(deepcc_model, new_fs_single)
```
Note: You can generate customized reference expression profile from your previous data or public data, which is the same(similar) cancer type and platform. Alternatively, you can use pre-defined reference in DeepCC by passing the cancer type (in the format of TCGA cancer types)

### Additional tools
- `cross_validattion` performs cross validation.
- `get_gene_sets` get gene sets list from MSigDB.
- `vis_sample` visualize samples
```
#use the deep feature for visualization of classification. 
#Assume deep feature stored in df, labels gotten form `get_DeepCC_label` stored in labels, 
#and color defined by yourself which should have the same length with the levels of labels.
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
vis_samples(df, labels, color)
```
## Pre-defined data in DeepCC
### Colorectal cancer model trained with public data sets 
Provide Colorectal cancer DeepCC model used in DeepCC_online. One is trainded with TCGA-COADREAD data set with 626 samples. Another is trainded with CRCSC data sets. (**Version: 20200901_v0.1.1**)

#### baidunetdisk
Link: [https://pan.baidu.com/s/10ogQPsNa4SksJadT8_ceHA](https://pan.baidu.com/s/10ogQPsNa4SksJadT8_ceHA)   &nbsp;&nbsp;&nbsp;&nbsp; Password: ugeg

To load CRC DeepCC model in R:
```
load_DeepCC_model <- function(prefix){
  load(file = paste0(prefix, ".RData"))
  classifer <- keras::load_model_hdf5(filepath =paste0(prefix, ".hdf5"))
  list(classifier = classifer, levels = levels)
}

# for CRC model trained with TCGA data
CRC_TCGA <- load_DeepCC_model("CRC_TCGA")
```
### DeepCC_online models
All DeepCC models used in DeepCC_online are available in [github](https://github.com/zero19970/deepcc_model):
```
# use git in terminal
git clone git@github.com:zero19970/deepcc_model.git
```
### List of functional gene DeepCC
By default, DeepCC will use MSigDB v7 (22, 596 gene sets) to generate functional spectra. You can also use MSigDB v5.0 (10, 348 gene sets) and MSigDB v7.0 (17, 779 gene sets).
```
# get MSigDBr from R package `msigdbr`
library(msigdbr)
MSigDBr <- get_msigdbr()
```

### Pre-defined reference
In DeepCC we prepared average expression profiles of each cancer types in TCGA project as reference. To use them, just use the TCGA identifier (COADREAD, BRCA, OV, etc.) to indicate the cancer type.

Note: if your single sample is microarray data, we strongly suggest you turn the parameter `inverseRescale` on, since TCGA is RNA-Seq data, compared with microarray. `inverseRescale` can overcome this distribution differently a little bit, but it's better to generate your customized reference.
