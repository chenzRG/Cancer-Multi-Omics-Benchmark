\name{plotGS}
\alias{plotGS}

\title{
  Plot the gene set space
}
\description{
  Plot the gene set space of objects of "moa" and "mgsa"
}
\usage{
  plotGS(x, axes=1:2, center.only=FALSE, topN=1, data.pch=20, data.col=1, highlight.col = 2,
         label=NULL, label.cex=1, layout=NULL, ...)
}

\arguments{
  \item{x}{
    An object of class \code{\link{mgsa-class}} or \code{\link{moa.sup-class}}
}
  \item{axes}{
    An integer vector in the length 2 to indicate the axes to be drawn.
}  
  \item{center.only}{
    A logical to indicate whether the separate gene set spaces from each of the data set should be plotted. Default is FALSE.
}
  \item{topN}{
    An integer specify N gene set from the most positive and negative end of axes to be labeled
}
  \item{data.pch}{
    The shape for plotting each data set. This argument is passed to points function, so only used when separate gene set spaces are plotted (i.e. center.only = FALSE).
}
  \item{data.col}{
    The col for plotting each data set. This argument is passed to points function, so only used when separate gene set spaces are plotted (i.e. center.only = FALSE).
}
  \item{highlight.col}{
    The color used to highlight the selected gene sets
}
  \item{label}{
    Either a character vector or NULL (default). The character vector should be the name of some gene sets want ot be labeled.
}
  \item{label.cex}{
    Passed to \code{\link{text}} function to adjust the the labels 
}
  \item{layout}{
    A matrix passed to the \code{layout} function.
}
  \item{\dots}{
    Other arguments passed to \code{\link{points}}
}
}

\details{
  This is a convenience function to explore the gene set space so not very flexible. For customized plot, please use the object of \code{data@coord.comb} and \code{data@coord.sep}.
}

\value{
    If assign to variable, A \code{list} of selected/highlighted gene set at 
    the (positve and negative) end of each axis will be returned.
}

\author{
  Chen Meng
}

\examples{
  
  # library(mogsa)
  # loading gene expression data and supplementary data
  data(NCI60_4array_supdata)
  data(NCI60_4arrays)
  mgsa <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=9,
                proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)

  plotGS(mgsa, center.only = TRUE, topN=5)
  res <- plotGS(mgsa, center.only = FALSE, data.pch=1:4, data.col=1:4)
  res
}
