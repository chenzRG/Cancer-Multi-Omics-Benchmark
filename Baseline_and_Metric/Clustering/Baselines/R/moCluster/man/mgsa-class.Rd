\name{mgsa-class}
\Rdversion{1.1}
\docType{class}
\alias{mgsa-class}

\title{Class \code{"mgsa"}}
\description{
	mgsa class here.
}
\section{Objects from the Class}{
	Objects can be created by calls of the form \code{new("mgsa", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{ call }
    \item{\code{moa}:}{Object of class \code{moa} }
    \item{\code{sup}:}{Object of class \code{moa.sup} }
  }
}
\section{Methods}{
	\describe{
  		\item{combine}{\code{signature(x = "mgsa", y = "mgsa")}}{
  			To combine two objects of class "mgsa"
  		}
  		This function could only be used to combine two "mgsa" objects, using "Reduce" function to combine more.

      \item{show}{\code{signature(x = "moa", y = "missing")}:
        show the "mgsa" result.
      }
  	}
}

\author{
	Chen Meng
}

\seealso{
	\code{\linkS4class{moa}} and \code{\linkS4class{moa.sup}}
}
\examples{
  showClass("mgsa")
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
  dim(getmgsa(mgsa1, "fac.scr"))
  dim(getmgsa(mgsa2, "fac.scr"))
  dim(getmgsa(mgsa_comb, "fac.scr"))
}
\keyword{classes}
