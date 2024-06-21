\name{toMoa}
\alias{toMoa}

\title{
convert \code{\link{mbpca}} result to \code{\link{moa-class}}
}
\description{
An internal function called by \code{\link{mbpca}}.
}
\usage{
toMoa(data, x, call)
}

\arguments{
  \item{data}{
  The preprocessed data in \code{\link{mbpca}}
}
  \item{x}{
  The object calculated in \code{\link{mbpca}}
}
  \item{call}{
  The call of mbpca
}
}

\value{
  An object of \code{moa-class}.
}

\author{
  Chen Meng
}
