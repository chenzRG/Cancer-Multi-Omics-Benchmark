\name{combine-methods}
\docType{methods}
\alias{combine}
\alias{combine-methods}
\alias{combine,mgsa,mgsa-method}

\title{
  Combine two objects of class \code{mgsa} into one.
}

\usage{
  combine(x, y, ...)
}

\arguments{
  \item{x}{ one mgsa object  }
  \item{y}{ another mgsa object }
  \item{...}{ ignored. Only two mgsa objects could be combined, using "Reduce"
  to combine more than two sets. }
}

\description{
  This function could only be used to combine two "mgsa" objects at present; using "Reduce" function to combine more.
}

\section{Methods}{
  \describe{
    \item{\code{signature(x = "mgsa", y = "mgsa")}}{ 
      To combine two objects of \code{mgsa}.
    }
    This function could only be used to combine two "mgsa" objects; using "Reduce" function to combine more.
  }
}

\value{
    A combined object of class \code{mgsa} will be returned. 
}

\examples{
  # library(mogsa)
  # loading gene expression data and supplementary data
  data(NCI60_4array_supdata)
  data(NCI60_4arrays)
  # split gene set annotation into two sets.
  sup1 <- lapply(NCI60_4array_supdata, function(x) x[, 1:10])
  sup2 <- lapply(NCI60_4array_supdata, function(x) x[, -(1:10)])
  # project two sets of annotation
  mgsa1 <- mogsa(x = NCI60_4arrays, sup=sup1, nf=9,
                proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
  mgsa2 <- mogsa(x = NCI60_4arrays, sup=sup2, nf=9,
                 proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
  # combine two indenpendent mgsa sets
  mgsa_comb <- combine(mgsa1, mgsa2)
  dim(getmgsa(mgsa1, "score"))
  dim(getmgsa(mgsa2, "score"))
  dim(getmgsa(mgsa_comb, "score"))
}

\keyword{combine}
\keyword{mogsa}
\keyword{mgsa-class}
