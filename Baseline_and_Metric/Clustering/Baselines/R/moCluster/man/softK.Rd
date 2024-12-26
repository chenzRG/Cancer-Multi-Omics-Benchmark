\name{softK}
\alias{softK}

\title{
  Soft-thresholding operator
}
\description{
  Weighted soft-thresholding operator, which is called by \code{\link{mbpca}}.
}
\usage{
softK(x, k, w = 1, pos = FALSE)
}

\arguments{
  \item{x}{
  A numerical vector
}
  \item{k}{
  Number of non-zero elements want to keep
}
  \item{w}{
  weight for each element. The actual thresholding is base on x*w, the default setting equals to ordinary soft thresholding. 
}
  \item{pos}{
  A logical value, if only positive values are retained.
}
}

\value{
    A thresholded numerical vector
}

\author{
    Chen Meng
}

\examples{
v <- rnorm(10)
softK(v, k = 2)
}

\keyword{soft threshold}
\keyword{soft thresholding}
\keyword{weighted soft threshold}
\keyword{weighted soft thresholding}