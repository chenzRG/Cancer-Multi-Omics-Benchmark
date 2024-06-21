\name{wsvd}
\alias{wsvd}

\title{
    Weighted singular value decomposition (SVD)
}
\description{
    The weighted version of singular value decomposition.
}
\usage{
    wsvd(X, D1 = diag(1, nrow(X)), D2 = diag(1, ncol(X)))
}

\arguments{
  \item{X}{
    A numeric matrix whose wSVD decomposition is to be computed.
}
  \item{D1}{
    A square matrix or vector. The left constraint/weight matrix (symmetric and positive in diagonal). The dimension of D1 should be the same with the number of rows in X. A vector input will be converted to a diagnal matrix.
}
  \item{D2}{
    A square matrix or vector. The right constraint/weight matrix (symmetric, positive in diagonal). The dimension of D1 should be the same with the number of columns in X. A vector input will be converted to a diagnal matrix.
}
}
\details{
    The weighted version of generalized singular value decomposition (SVD) of matrix 
    A = UDV' 
    with the constraints 
    U'D1U = I and V'D2V = I
    D1 and D2 are two matrices express constraints imposed on the rows and the columns of matrix A.
}
\value{
    d - singular values

    u - left singular vectors

    v - right singular vectors

    D1 - the left weight matrix (directly from input)
    
    D2 - the right weight matrix (directly from input)
}
\references{
    Herve Abdi. Singular Value Decomposition (SVD) and Generalized Singular Value Decomposition (GSVD)
    http://www.utdallas.edu/~herve/Abdi-SVD2007-pretty.pdf
}
\author{
    Chen Meng
}

\seealso{
 svd
}

\examples{
    set.seed(56)
    m <- matrix(rnorm(15), 5, 3)
    wl <- rnorm(5)
    wr <- runif(3)
    s <- wsvd(X=m, D1=wl, D2=wr)
    # t(s$u) %*% diag(wl) %*% s$u
    # t(s$v) %*% diag(wr) %*% s$v
    # all.equal(m, as.matrix(s$u) %*% diag(s$d) %*% t(s$v))
}

\keyword{ SVD }
\keyword{ generalized SVD }
\keyword{ weighted SVD }
