\name{processOpt}
\alias{processOpt}

\title{
preprocessing of input data in \code{\link{mbpca}}.
}
\description{
An internal function called by \code{\link{mbpca}}.
}
\usage{
processOpt(x, center = TRUE, scale = FALSE, option = c("lambda1", "inertia", "uniform"))
}

\arguments{
  \item{x}{
  A list of matrices, rows are observations and columns are variables
}
  \item{center}{
  A logical variable indicates whether columns should be centered
}
  \item{scale}{
  A logical variable indicates whether columns should be scaled
}
  \item{option}{
  A charater string could be one of c("lambda1", "inertia", "uniform") to indicate
  how the different matrices should be normalized. If "lambda1", the matrix is 
  divided by its the first singular value, if "inertia", the matrix is divided by
  its total inertia (sum of square), if "uniform", none of them would be done.
}
}

\value{
  A \code{list} of normalized matrix.
}

\author{
  Chen Meng
}

