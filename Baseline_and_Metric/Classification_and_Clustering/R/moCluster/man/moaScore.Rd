\name{moaScore}
\alias{moaScore}

\title{
  Extract global scores from an object of class \code{\link{moa-class}}.
}
\description{
  Extract global scores from an object of class \code{\link{moa-class}}.
}
\usage{
moaScore(moa)
}

\arguments{
  \item{moa}{
	An object of class \code{\link{moa-class}}
}
}

\value{
	A matrix of global score
}

\author{
	Chen Meng
}


\seealso{
\code{\link{moaCoef}}
}

\examples{
# see examples in \code{\link{mbpca}}

data("NCI60_4arrays")
moa <- mbpca(NCI60_4arrays, ncomp = 10, k = "all", method = "globalScore", option = "lambda1", 
             center=TRUE, scale=FALSE)

genes <- moaCoef(moa)
scr <- moaScore(moa)

}
