\name{distMoa}
\alias{distMoa}

\title{
  Calculate the distance matrix from an object of class \code{\link{moa-class}}.
}
\description{
  A convenient function to calculate the distance matrix from 
  an object of class \code{\link{moa-class}}.
}
\usage{
distMoa(x, nf = NA, tol = 1e-05, method = "euclidean", 
  diag = FALSE, upper = FALSE, p = 2)
}

\arguments{
  \item{x}{
  An object of class \code{\link{moa-class}}.
}
  \item{nf}{
  Integer; the number of component used to calculate the distance. Default setting (NA) will keep all the axes.
}
  \item{tol}{
  Numerical; the tolerance of component with low variance.
}
  \item{method}{
  passed to function \code{dist}
}
  \item{diag}{
  passed to function \code{dist}
}
  \item{upper}{
  passed to function \code{dist}
}
  \item{p}{
  passed to function \code{dist}
}
}

\value{
An object of class \code{dist}, see function "dist".
}

\author{
Chen Meng
}

\examples{
# see examples in \code{\link{mbpca}}

data("NCI60_4arrays")
moa <- mbpca(NCI60_4arrays, ncomp = 10, k = "all", method = "globalScore", option = "lambda1", 
             center=TRUE, scale=FALSE)

dst <- distMoa(moa)


}

