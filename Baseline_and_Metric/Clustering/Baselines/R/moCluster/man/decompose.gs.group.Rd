\name{decompose.gs.group}
\alias{decompose.gs.group}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Data-wise or PC-wise decomposition of gene set scores for all observations.
}
\description{
  Data-wise or PC-wise decomposition of gene set scores (GSS) across all observations. The predefined group/cluster information should be given so that the mean decomposed GSSs for each group are returned and plotted.
}
\usage{
  decompose.gs.group(x, gs, group, decomp = "data", nf = 2, x.legend = "bottomleft", 
    y.legend = NULL, plot = TRUE, main = NULL, ...)
}

\arguments{
  \item{x}{
    An object of class \code{\link{mgsa-class}} or \code{\link{moa.sup-class}}
}
  \item{gs}{
    The gene set want to exam.
}
  \item{group}{
    An vector or factor to indicate the group of observations, such as clusters. 
    See examples.
}
  \item{decomp}{
    A charater string either "data" or "pc" to indicate how the gene set scores
    should be decomposed (with respect to data or PC.
}
  \item{nf}{
    The number of axes/PCs to be calculated and plotted.
}
  \item{x.legend}{
    Used to control the position of legends.
}
  \item{y.legend}{
    Used to control the position of legends.
}
  \item{plot}{
    A logical indicates if a plot should be drawn.
}
  \item{main}{
    The main title of plot.
}
  \item{\dots}{
    Other arguments passed to \code{barplot}.
}
}

\details{
  This function could be used when the number of observation is large and there
  are cluster/group information is available. In this case, the means of 
  decomposed gene set scores over each group is calculated. The vertical bar
  on the end of each bar indicates the 95\% confident interval of the means.
}

\value{
  Return nothing or a matrix depends on how argument \code{plot} is set.
}

\references{
  TBA
}
\author{
  Chen Meng
}

\seealso{
  See Also \code{\link{decompose.gs.ind}}
}

\examples{
  # library(mogsa)
  # loading gene expression data and supplementary data
  data(NCI60_4array_supdata)
  data(NCI60_4arrays)

  # using a list of data.frame as input
  mgsa <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=9,
                proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)

  colcode <- as.factor(sapply(strsplit(colnames(NCI60_4arrays$agilent), split="\\\."), "[", 1))
  decompose.gs.group(x = mgsa, gs = 2, group = colcode, decomp = "data", plot = TRUE)
  decompose.gs.group(x = mgsa, gs = 2, group = colcode, decomp = "pc", nf = 3, plot = TRUE)
}

