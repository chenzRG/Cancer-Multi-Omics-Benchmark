\name{moaCoef}
\alias{moaCoef}

\title{
Extract the loadings/coefficients from an object of class \code{\link{moa-class}}.
}
\description{
Extract the loadings/coefficients from an object of class \code{\link{moa-class}}.
}
\usage{
moaCoef(moa)
}

\arguments{
  \item{moa}{
  An object of class \code{\link{moa-class}}.
}
}

\value{
    It returns a list consist of two components:

    \code{coefMat} - the loading matrix

    \code{nonZeroCoef} - it is a \code{list} of \code{data.frame} to list the non-zero coefficient
    variable in each of loading vectors and data sets. The element names are in a format as

    "xxxx.yy.zzz"

    xxxx - are the data names, tells the data set where a varirable is from

    yy - the number of Axes, for example, "V1" indicate the variable has a non-zero coefficient
    in the first loading vector.

    zzz - could be either "pos" (coefficient >0) or "neg" (coefficient < 0)

    The \code{data.frame} has two columns, the first column is the ID of a variable the second 
    column is the coefficient/loading.
}

\author{
    Chen Meng
}

\seealso{
\code{\link{moaScore}}
}

\examples{
# see examples in \code{\link{mbpca}}

data("NCI60_4arrays")
moa <- mbpca(NCI60_4arrays, ncomp = 10, k = "all", method = "globalScore", option = "lambda1", 
             center=TRUE, scale=FALSE)

genes <- moaCoef(moa)
scr <- moaScore(moa)

}

