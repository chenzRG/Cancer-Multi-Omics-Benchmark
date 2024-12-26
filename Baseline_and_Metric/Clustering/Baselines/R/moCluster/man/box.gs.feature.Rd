\name{box.gs.feature}
\alias{box.gs.feature}

\title{
  boxplot of gene set variables across all samples.
}
\description{
  boxplot to show the variables (e.g. gene expression) of a gene set across all samples. 
}
\usage{
  box.gs.feature(x, gs, moa = NULL, col = 1, layout = NULL, plot = TRUE, obs.order = NULL, ...)
}

\arguments{
  \item{x}{
    An object of calss \code{\link{mgsa-class}} or \code{\link{moa.sup-class}}
  }
  \item{gs}{
    Gene set want to be explored
  }
  \item{moa}{
    An obejct of class \code{\link{moa}}. It is required if x is an object of 
    class \code{\link{moa.sup-class}}
  }
  \item{col}{
    The coler code for samples
  }
  \item{layout}{
    The layout control, see examples.
  }
  \item{plot}{
    A logical indicates whether the result should be ploted. If FALSE, a list of expression 
    matrix of the gene set genes is returned. Otherwise nothing returned.
  }
  \item{obs.order}{
    Can be used to reorder the martrix, could be used when clustering result is available.
  }
  \item{\dots}{
    The arguments passed to \code{\link{boxplot}}
  }
}
\details{
  This is a convenient function used to explore the expression of a set of features/genes
}
\value{
  Do not return anything (plot=TRUE) or return a list of matrix (plot=FALSE) depends on 
  plot arugment.
}

\author{
  Chen meng
}

\examples{
  # library(mogsa)
  # loading gene expression data and supplementary data
  data(NCI60_4array_supdata)
  data(NCI60_4arrays)
  mgsa <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=9,
                proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)

  allgs <- colnames(NCI60_4array_supdata[[1]])
  colcode <- as.factor(sapply(strsplit(colnames(NCI60_4arrays$agilent), split="\\\."), "[", 1))
  a <- box.gs.feature(x=mgsa, gs=allgs[5], type=3, col=colcode, plot=FALSE)
  box.gs.feature(x=mgsa, gs=allgs[5], type=3, col=colcode, plot=TRUE, layout=matrix(1:4, 2, 2))
}
