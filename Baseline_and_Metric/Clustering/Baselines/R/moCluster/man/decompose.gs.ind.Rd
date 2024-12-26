\name{decompose.gs.ind}
\alias{decompose.gs.ind}

\title{
  Data-wise or PC-wise decomposition of gene set scores for a single observation.
}

\description{
  Barplot of decomposed gene set scores, either with respect to datasets
  or axes.
}

\usage{
  decompose.gs.ind(x, gs, obs, type = 3, nf = 2, plot=TRUE, col.data = NULL, 
    col.pc = NULL, legend = TRUE)
}

\arguments{
  \item{x}{
    An object of class \code{\link{mgsa-class}} or \code{\link{moa.sup-class}}
  }
  \item{gs}{
    The gene set want to exam.
  }
  \item{obs}{
    The observations want to exam.
  }
  \item{type}{
    Which type of plot.
    type=1 - the data-pc mode; 
    type=2 - the pc-data mode; 
    type=3 - both. See detail.
  }
  \item{nf}{
    The number of axes/PCs to be calculated and plotted.
  }
  \item{plot}{
    A logical indicates if a plot should be drawn
  }
  \item{col.data}{
    The bar color of datasets
  }
  \item{col.pc}{
    The bar color of PCs
  }
  \item{legend}{
    A logical if legend should be shown
  }
}

\details{
  type=1 (the data-pc mode), the axes/PCs are represented as the narrow bars with different
  colors and the background wide bars behind narrow bars are gene set scores for datasets, 
  which is calculated from the sum of all underlying individual axes/PC scores.
  When type=2 (the pc-data mode) the interpreation of narrow and wide bars are in the other way around. If type=3, both are shown.

  This function could only be used to check the decomposition of gene set scores 
  of a single observation. So the function is not efficent when the number of 
  observation is large. Another function \code{\link{decompose.gs.group}}, could be 
  used in this case, particularly when the cluster information of the observation 
  panel is available.
}

\value{
  Return nothing or a matrix depends on how argument \code{plot} is set.
}

\references{
  TBA
}

\author{
  Chen Meng
}

\seealso{
  See Also as \code{\link{decompose.gs.group}}
}

\examples{
  # library(mogsa)
  # loading gene expression data and supplementary data
  data(NCI60_4array_supdata)
  data(NCI60_4arrays)
  mgsa <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=9,
                proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)

  allgs <- colnames(NCI60_4array_supdata[[1]])
  # plot
  decompose.gs.ind(x=mgsa, gs=allgs[5], obs="BR.MDA_MB_231", type=2, nf=5)
  # or
  decompose.gs.ind(x=getmgsa(mgsa, "sup"), gs=allgs[5], obs="BR.MDA_MB_231", type=3, nf=5)
}


