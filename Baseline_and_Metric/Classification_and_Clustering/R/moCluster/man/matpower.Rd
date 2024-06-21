\name{matpower}
\alias{matpower}

\title{
 compute the power of a matrix
}
\description{
 the power of a matrix
}
\usage{
matpower(x, n, nf = min(dim(x)), tol = 1e-07)
}

\arguments{
  \item{x}{
  	a numerical matrix object that the power of which should be calculated
}
  \item{n}{
  	The matrix to the power of
}
  \item{nf}{
	The number of axes kept in the calculation of SVD and reconstruction
}
  \item{tol}{
	The tolerance of the axis, singular vectors with singular value lower than tol will be ignored in the reconstruction.
}
}
\details{
	The power of a matrix is calculated in two steps:
	decompostion step: x=UDV'
	and the reconstruction step: x^n=U*D^n*V'
	In the reconstruction, the singular vectors with a singular value more than
	tol are kept.
}
\value{
	A matrix x^n
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
	Chen Meng
}
\note{
	Called by the wsvd function.
}

\seealso{
	See Also \code{\link{wsvd}}
}
\examples{
  set.seed(56)
  m <- matrix(rnorm(15), 5, 3)
  s <- matpower(m, 2)
  s <- matpower(m, -2)
}

