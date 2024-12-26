\name{nipalsSoftK}
\alias{nipalsSoftK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    NIPALS algorithm with soft thresholding operator
}
\description{
    An internal function called by \code{\link{mbpca}}.
}
\usage{
nipalsSoftK(x, maxiter, k)
}

\arguments{
  \item{x}{
  The input matrix, rows are observations, columns are variables
}
  \item{maxiter}{
  # of maximum interation the algorithm can run
}
  \item{k}{
  The number (>=1) or proportion (<1) of variables want to keep. 
  It could be a single value or a vector has the same length as x so 
  the sparsity of individual matrix could be different.
}
}

\value{
    an \code{list} object contains the following elements:

    \code{tb} - the block scores 

    \code{pb} - the block loadings

    \code{t} - the global scores

    \code{w} - the wegihts of block scores to construct the global score.
}

\author{
    Chen Meng
}

\seealso{
\code{\link{msvd}}
}