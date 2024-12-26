\name{deflat}
\alias{deflat}

\title{
  deflat function used by \code{\link{mbpca}}
}
\description{
  An internal function called by \code{\link{mbpca}}.
}
\usage{
deflat(x, t, tb, pb, method = "globalScore")
}

\arguments{
  \item{x}{
  A \code{list} of \code{matrix} want to deflat
}
  \item{t}{
  The global scores returned by \code{\link{msvd}} or \code{\link{nipalsSoftK}}
}
  \item{tb}{
  The block scores returned by \code{\link{msvd}} or \code{\link{nipalsSoftK}}
}
  \item{pb}{
  The block loadings returned by \code{\link{msvd}} or \code{\link{nipalsSoftK}}
}
  \item{method}{
  A charater to specify the deflation strateg, could be one of 
  c("globalScore", "blockLoading", "blockScore").
}
}

\value{
A \code{list} of deflated \code{matrix}
}

\author{
Chen Meng
}
