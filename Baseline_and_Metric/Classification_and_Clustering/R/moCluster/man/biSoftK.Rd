\name{biSoftK}
\alias{biSoftK}

\title{
    NIPALS algorithm with soft thresholding operator on rows and columns
}

\description{
    An internal function called by \code{\link{mbpca}}.
}
\usage{
biSoftK(x, maxiter, kp, kt, weight.p, weight.t, pos = FALSE, unit.pb = TRUE, unit.tb = FALSE)
}

\arguments{
  \item{x}{
  The input matrix, rows are observations, columns are variables
}
  \item{maxiter}{
  Number of maximum interation the algorithm can run
}
  \item{kp}{
  The number (>=1) or proportion (<1) of variables want to keep. 
  It could be a single value or a vector has the same length as x so 
  the sparsity of individual matrix could be different.
} 
  \item{kt}{
  The number (>=1) or proportion (<1) of non-zero scores for obvservations. 
}
  \item{weight.p}{
  The weight of variables. It could be 1) a vector has the same length as x, one value for each table/block; 
  2) one number, all variables share the same weight or 3) a list of vectors, the length of each vector should
  be the same with the columns numbers of the corresponding table/block, so every variables has a unique weight.
}
  \item{weight.t}{
  The weight for observation. For accepted values or formats, see weight.p.
}
  \item{pos}{
  Logical value, if only non-negaitve values in the loading and score vectors.
}
  \item{unit.pb}{
  Logical value, whether the length of table/block loading should be unit length. 
}
  \item{unit.tb}{
  Logical value, whether the length of table/block score should be unit length. 
}
}

\details{
  This function also use the NIPALS algorithm, but it generalized nipalsSoftK from several aspects:
  1. Allowing sparsity on both columns and rows of matrices
  2. Allowing weights for columns and rows
  3. Allowing loading and/or score vectors of blocks to be unit length
  4. Allowing only positive number in loading and score vectors
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