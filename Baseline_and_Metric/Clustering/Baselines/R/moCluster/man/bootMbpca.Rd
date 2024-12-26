\name{bootMbpca}
\alias{bootMbpca}

\title{
Bootstrap mbpca to estimate the coherence of different data sets
}
\description{
Bootstrap mbpca to estimate the coherence of different data sets and estimate the 
number of components should be included in an analysis.
}
\usage{
bootMbpca(moa, mc.cores = 1, B = 100, replace = TRUE, 
  resample = c("sample", "gene", "total"), log = "y", ncomp = NULL, method = NULL, 
  maxiter = 1000, svd.solver = c("svd", "fast.svd", "propack"), plot = TRUE)
}

\arguments{
  \item{moa}{
  An object of \code{\link{moa}} returned by \code{\link{mbpca}}.
}
  \item{mc.cores}{
  Integer; number of cores used in bootstrap. This value is passed to function \code{mclapply}
}
  \item{B}{
  Integer; number of bootstrap
}
  \item{replace}{
  Logical; sampling with or without replacement
}
  \item{resample}{
  Could be one of "sample", "gene" or "total". "sample" and "gene" means sample-wise 
  and variable-wise resampling, repectively. "total" means total resampling.
}
  \item{log}{
  Could be "x", "y" or "xy" for plot log axis
}
  \item{ncomp}{
  Passed to function \code{\link{mbpca}}. In most of cases, user do not need to specify this argument
  because it could be inferred from \code{moa}.
}
  \item{method}{
  Passed to function \code{\link{mbpca}}.In most of cases, user do not need to specify this argument
  because it could be inferred from \code{moa}.
}
  \item{maxiter}{
  Passed to function \code{\link{mbpca}}.In most of cases, user do not need to specify this argument
  because it could be inferred from \code{moa}.
}
  \item{svd.solver}{
  Passed to function \code{\link{mbpca}}.In most of cases, user do not need to specify this argument
  because it could be inferred from \code{moa}.
}
  \item{plot}{
  Logical; whether the result should be plotted.
}
}
\details{
  Bootstrap method were used to determine the components that are presenting significant concordant structure between datasets.
}
\value{
It returns a matrix, columns are eigenvalues for different components. Each rows is a bootstramp sample.
}

\author{
Chen Meng
}

\examples{
# see examples in \code{\link{mbpca}}
}

