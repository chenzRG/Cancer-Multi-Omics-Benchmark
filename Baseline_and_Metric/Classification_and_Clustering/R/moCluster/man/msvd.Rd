\name{msvd}
\alias{msvd}

\title{
    SVD based algorithm to calculate block Score and global scores for \code{\link{mbpca}}.
}
\description{
    An internal function called by \code{\link{mbpca}}. It returns the result comparable 
    with nipalsSoftK, but way faster since it uses the SVD algorithm. No sparse opertors
    in this function.
}
\usage{
msvd(x, svd.sol = svd)
}

\arguments{
  \item{x}{
  The input matrix, rows are observations, columns are variables
}
  \item{svd.sol}{
  A function object to specify the preferred SVD solver, default is \code{svd}.
}
}

\value{
    an \code{list} object contains the following elements:

    \code{tb} - the block scores 

    \code{pb} - the block loadings

    \code{t} - the global scores

    \code{w} - the wegihts of block scores to construct the global scor
}

\author{
    Chen Meng
}

\seealso{
\code{\link{nipalsSoftK}}
}
