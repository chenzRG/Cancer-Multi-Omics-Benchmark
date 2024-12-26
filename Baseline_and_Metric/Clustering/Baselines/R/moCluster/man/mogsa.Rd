\name{mogsa}
\alias{mogsa}

\title{
    multiple omics data integration and gene set analysis
}

\description{
    The main function called by users, omics data analysis and gene set annotation. 
    A wrapper function of \code{\link{moa}} and \code{\link{sup.moa}}.
}

\usage{
  mogsa(x, sup, nf=NULL, factors = NULL, proc.row=NULL, w.data=NULL, w.row=NULL, statis=FALSE, ks.stat=FALSE, ks.B = 1000, ks.cores = NULL, p.adjust.method = "fdr")
}

\arguments{
    \item{x}{
      An object of class \code{list} or \code{\link{moa-class}}. A list would be a list of data frame. 
    }
    \item{sup}{
      An object of class \code{list} or \code{\link{moa.sup-class}}. A list would be a list of supplementary data.
    }
    \item{nf}{
      The number of principal components used to reconstruct, only used when x is a an object of \code{list}.
    }
    \item{factors}{
    The index of principal components used in the projection, used when non-consecutive PC to be included in the analysis.
    }
    \item{proc.row}{
      Preprocessing of rows. If x is a object of \code{list}, it is passed \code{moa}
    }
    \item{w.data}{
      Weights of datasets. If x is a object of \code{list}, it is passed \code{moa}
    }
    \item{w.row}{
      Weight of row. If x is a object of \code{list}, it is passed \code{moa}
    }
    \item{statis}{
      A logical indicates if statis algrithm should be used. If x is a object of \code{list}, it is passed \code{moa}
    }
    \item{ks.stat}{
      The logical indicates if the p-value should be calculated using K-S statistic (the method used in "ssgsea" in GSVA package). 
      Default is FALSE, which means using the z-score method. See \code{\link{sup.moa}}.
    }
    \item{ks.B}{
      An integer to indicate the number of bootstrapping samples to calculated the p-value of KS statistic.
    }
    \item{ks.cores}{
      An integer indicate the number of cores to be used in bootstrapping. It is passed to function \code{mclapply} in the 
      \code{parallel} package.
    }
    \item{p.adjust.method}{
    The method of p value adjustment, passed to \code{p.adjust} function.
    }
}

\details{
    A wrapper function of \code{\link{moa}} and \code{\link{sup.moa}}. 
}

\value{
    An object of class \code{\link{mgsa-class}}. 
}

\references{
    Preprint: Meng, C., Kuster, B., Peters, B., Culhane, AC., Moghaddas Gholami, A., moGSA: integrative single sample gene-set analysis of multiple omics data. doi: http://dx.doi.org/10.1101/046904
    Haenzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set variation analysis for microarray and RNA-Seq data. BMC Bioinformatics, 14:7, 2013.
    Barbie, D.A. et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature, 462(5):108-112, 2009.
}

\author{
    Chen Meng 
}

\note{
  This function will be changed to a generic function for "S4-style" programming.
}

\seealso{
    \code{\link{moa}} and \code{\link{sup.moa}}
}

\examples{
  # library(mogsa)
  # loading gene expression data and supplementary data
  data(NCI60_4array_supdata)
  data(NCI60_4arrays)

  # using a list of data.frame as input
  mgsa1 <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=9,
                 proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
  mgsa1x <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, factors = c(1,3,6),
                 proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
  # using moa as input
  ana <- moa(NCI60_4arrays, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
  smoa <- sup.moa(ana, sup=NCI60_4array_supdata, nf=3)
  mgsa2 <- mogsa(x = ana, sup=NCI60_4array_supdata, nf=9)
  mgsa3 <- mogsa(x = ana, sup=smoa)
}

\keyword{ MVA }
\keyword{ supplementary data projection }
