\name{sup.moa}
\alias{sup.moa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Projecting supplementary tables on object of class \code{moa-class}.
}
\description{
    Projecting supplementary tables on \code{\link{moa-class}}
}
\usage{
    sup.moa(X, sup, nf = 2, factors = NULL, ks.stat=FALSE, ks.B = 1000, ks.cores = NULL, p.adjust.method = "fdr")
}

\arguments{
  \item{X}{
    An object of class \code{\link{moa-class}}
  }
  \item{sup}{
    A list of data.frames contains supplementary data.
  }
  \item{nf}{
    The number of principal components used in the projection.
  }
  \item{factors}{
    The index of principal components used in the projection, used when non-consecutive PC to be included in the analysis.
  }
  \item{ks.stat}{
    The logical indicates if the p-value should be calculated using K-S statistic (the method used in "ssgsea" in GSVA package). 
    Default is FALSE, which means using the z-score method. 
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
    Projecting supplementary tables on \code{\link{moa-class}}, for details see reference.
}

\value{
    An object of class \code{\link{moa.sup-class}}. 
}

\references{
    Herve Abdi, Lynne J. Williams, Domininique Valentin and Mohammed Bennani-Dosse. STATIS and DISTATIS: optimum multitable principal component analysis and three way metric multidimensional scaling. WIREs Comput Stat 2012. Volume 4, Issue 2, pages 124-167
    Haenzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set variation analysis for microarray and RNA-Seq data. BMC Bioinformatics, 14:7, 2013.
    Barbie, D.A. et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature, 462(5):108-112, 2009.
}
\author{
    Chen Meng
}

\examples{
    # library(mogsa)
    # loading gene expression data and supplementary data
    data(NCI60_4array_supdata)
    data(NCI60_4arrays)
    # check the dimension of each supplementary data to see how many gene set annotated the data
    sapply(NCI60_4array_supdata, dim)
    # run analysis
    ana <- moa(NCI60_4arrays, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
    plot(ana, value="eig")
    # projectin supplementary data
    smoa <- sup.moa(ana, sup=NCI60_4array_supdata, nf=3)
    # heatmap visualize the gene set scores
    heatmap(slot(smoa, "score"))
}

\keyword{ data projection }
\keyword{ supplementary data}
