\name{annotate.gs}
\alias{annotate.gs}

\title{
  Summary annotation information of a gene set
}

\description{
  Retrive variables/features (genes) mapped to the annotated data sets in a gene set. Also returns the the information about presence and absence of a feature for a specific data set.
}

\usage{
  annotate.gs(mgsa, gs)
}

\arguments{
  \item{mgsa}{
    An object of class \code{\link{mgsa-class}} or \code{\link{moa.sup-class}}.
}
  \item{gs}{
    The name of a geneset
}
}

\value{
  This function returns a data.frame. 
  The first column shows the name of features.
  The last column is for the count of how many data sets has
  the corresponding features. Columns in the middle contains
  logical value indicating whether a feature is presented in a particular 
  data set.
}

\author{
  Chen Meng
}

\seealso{
 see \code{\link{GIS}}
}
\examples{
  # library(mogsa)
  # loading gene expression data and supplementary data
  data(NCI60_4array_supdata)
  data(NCI60_4arrays)
  mgsa <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=9,
                proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
  allgs <- colnames(NCI60_4array_supdata[[1]])
  annotate.gs(mgsa, allgs[1])
}

