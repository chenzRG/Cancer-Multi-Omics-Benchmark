Cancer Integration via Multikernel Learning (**CIMLR**)
=======================================================

**OVERVIEW**

In this repository we provide implementations in both R and Matlab of *CIMLR* (https://www.nature.com/articles/s41467-018-06921-8). This method was originally applied to multi-omic cancer data, but it is in principle capable of effectively and efficiently learning similarities in all the contexts where diverse and heterogeneous statistical characteristics of the data make the problem harder for standard approaches. 

The R branch of the repository provides the tool implemented in R, while the Matlab branch of the repository provides the Matlab implementation. Some example data are also provided, but those data are reduced versions of the original ones and should be used purely as examples and not considered as replacements of the ones provided in the respective publications. 

**CIMLR**

Outcomes for cancer patients vary greatly even within the same tumor type, and characterization of molecular subtypes of cancer holds important promise for improving prognosis and personalized treatment. This promise has motivated recent efforts to produce large amounts of multidimensional genomic ('multi-omic') data, but current algorithms still face challenges in the integrated analysis of such data. Here we present Cancer Integration via Multikernel Learning (CIMLR), a new cancer subtyping method that integrates multi-omic data to reveal molecular subtypes of cancer. We apply CIMLR to multi-omic data from 36 cancer types and show significant improvements in both computational efficiency and ability to extract biologically meaningful cancer subtypes. The discovered subtypes exhibit significant differences in patient survival for 27 of 36 cancer types. Our analysis reveals integrated patterns of gene expression, methylation, point mutations and copy number changes in multiple cancers and highlights patterns specifically associated with poor patient outcomes. 

**CITATION**

The latest version of the manuscript related to *CIMLR* is published on Nature Communications and can be found at https://www.nature.com/articles/s41467-018-06921-8. 

When using the tool, please cite: Ramazzotti, Daniele, et al. "Multi-omic tumor data reveal diversity of molecular mechanisms that correlate with survival." Nature Communications (2018); doi:10.1038/s41467-018-06921-8. 

**RUNNING CIMLR R IMPLEMENTATION**

The R version of *CIMLR* can be installed from Github. To do so, we need to install the R packages *CIMLR* depends on and the devtools package. 

First we run an R session and we execute the following commands. 

```r
# run this commands only if the following R packages are not already installed
install.packages("devtools", dependencies = TRUE)
install.packages("Matrix", dependencies = TRUE)
```

Now we can install and run *CIMLR* as follows: 

```r
# install CIMLR from Github
library("devtools")
install_github("danro9685/CIMLR", ref = 'R')

# load CIMLR library
library("CIMLR")
```

We provide a manual and an example of how to run CIMLR on multi-omic cancer data in the [doc folder](https://github.com/danro9685/CIMLR/tree/R/doc). 

**RUNNING CIMLR MATLAB IMPLEMENTATION**

The Matlab version of *CIMLR* is available in the Matlab branch of the repository. This version can be directly used through MATLAB by downloading the code in the repository. 

The scripts *Estimate_Number_of_Clusters_CIMLR.m* and *CIMLR.m* implement heuristics to estimate the best number of clusters from data and the code to run *CIMLR* respectively. The script *CIMLR_Feature_Ranking.m* provides code to assess which features were the most important for the clustering. 

A demo of the whole *CIMLR* analysis is provided in the script *main_example_lower_grade_gliomas.m* for a dataset of 282 patients affected by lower-grade glioma. The original data were published in Cancer Genome Atlas Research Network. "Comprehensive, integrative genomic analysis of diffuse lower-grade gliomas." New England Journal of Medicine 372.26 (2015): 2481-2498. 

**DEBUG**

Please feel free to contact us if you have problems running our tool at daniele.ramazzotti1@gmail.com or luca.desano@gmail.com. 
