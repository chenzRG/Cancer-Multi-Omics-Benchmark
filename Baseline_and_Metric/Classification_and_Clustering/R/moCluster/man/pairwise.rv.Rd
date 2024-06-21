\name{pairwise.rv}
\alias{pairwise.rv}

\title{
    pairwise RV coefficients.
}
\description{
    Calculating pairwise RV coefficients for a list of matrices or data.frame.
}
\usage{
    pairwise.rv(data.list, match="col")
}

\arguments{
  \item{data.list}{
    A list of data.frame or matrix, either rows or columns in each data set should be matched.
  }
    \item{match}{
    Whether columns or rows of data.frame/matrix should be matched.
  }
}
\details{
    The RV coefficient for each pair of matrices is calculated as 
        Rv = trace(XX'YY')/sqrt(trace(XX'XX')*trace(YY'YY'))
}

\note{
    The variable in matrices are not automatically centered or scaled in this function. So these step may need to be performed before calling this function.
}

\value{
    The function will return a matrix containing the pairwise RV coefficients.
}
\references{
    Robert, P.; Escoufier, Y. (1976). A Unifying Tool for Linear Multivariate Statistical 
    Methods: The RV-Coefficient. Applied Statistics 25 (3): 257-265.
}

\author{
    Chen Meng
}

\examples{
    data(NCI60_4arrays)
    pairwise.rv(NCI60_4arrays)
}

\keyword{ RV coefficent }

