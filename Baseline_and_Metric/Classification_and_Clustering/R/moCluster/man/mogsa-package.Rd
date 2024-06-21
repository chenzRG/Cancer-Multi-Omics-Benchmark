\name{mogsa-package}
\alias{mogsa-package}
% \alias{mogsa}
\docType{package}
\title{
Multiple omics clustering and gene set analysis
}
\description{
Modern "omics" technologies enable quantitative monitoring of the abundance of various biological molecules 
in a high-throughput manner, accumulating an unprecedented amount of quantitative information on a genomic 
scale. Gene set analysis is a particularly useful method in high throughput data analysis since it can 
summarize single gene level information into the biological informative gene set levels. This package
provide a method do the gene set analysis based on multiple omics data that describing the same set of 
observations/samples.
}
\details{
\tabular{ll}{
Package: \tab mogsa\cr
Type: \tab Package\cr
Version: \tab 1.3.1\cr
Date: \tab 2016-01-19\cr
License: \tab GPL-2\cr
Depends: \tab methods\cr
}
The main function in the package is "mogsa", see the function help manu for more details.
}
\author{
	Chen Meng
	Maintainer: Chen Meng <chen.meng@tum.de>
}
\references{
	Chen Meng, Dominic Helm, Martin Frejno, and Bernhard Kuster. 
  moCluster: Identifying Joint Patterns Across Multiple Omics Data
  Sets. Journal of Proteome Research 2016.


}

\examples{
  # library(mogsa)
  # loading gene expression data and supplementary data
  data(NCI60_4array_supdata)
  data(NCI60_4arrays)

  # using a list of data.frame as input
  mgsa1 <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=9,
                 proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
  # using moa as input
  ana <- moa(NCI60_4arrays, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
  smoa <- sup.moa(ana, sup=NCI60_4array_supdata, nf=3)
  mgsa2 <- mogsa(x = ana, sup=NCI60_4array_supdata, nf=9)
  mgsa3 <- mogsa(x = ana, sup=smoa)
}
