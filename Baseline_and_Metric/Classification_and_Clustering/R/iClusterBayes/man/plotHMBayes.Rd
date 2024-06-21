\name{plotHMBayes}
\alias{plotHMBayes}
\title{A function to generate heatmap panels sorted by integrated cluster assignment.}
\description{
  A function to generate heatmap panels sorted by integrated cluster assignment.
}
\usage{
plotHMBayes(fit, datasets, type = c("gaussian", "binomial", "poisson"),
    sample.order = NULL, row.order = NULL, sparse = NULL, 
    threshold = rep(0.5,length(datasets)), width = 5,
    scale = rep("none",length(datasets)), 
    col.scheme = rep(list(bluered(256)),length(datasets)),
    chr=NULL, plot.chr=NULL, cap=NULL)
}
\arguments{
  \item{fit}{A iClusterBayes object.}
  \item{datasets}{A list object of data matrices.}
  \item{type}{Types of data in the datasets.}
  \item{sample.order}{User supplied cluster assignment.}
  \item{row.order}{A vector of logical values each specify whether the genomic features in the corresponding data matrix should be reordered by similarity. Default is TRUE.}
  \item{sparse}{A vector of logical values each specify whether to plot the top cluster-discriminant features. Default is FALSE.}
  \item{threshold}{When sparse is TRUE, a vector of threshold values to include the genomic features on the heatmap.  Each data set should have a threshold.
  		  For each data set, a feature with posterior probability greater than the threshold will be included. Default value is 0.5 for each data set.}
  \item{width}{Width of the figure in inches}
  \item{scale}{A vector of logical values each specify whether data should be scaled. Default is FALSE.}
  \item{col.scheme}{Color scheme. Can use bluered(n) in gplots R package.}
  \item{chr}{A vector of chromosome number.}
  \item{plot.chr}{A vector of logical values each specify whether to annotate chromosome number on the left of the panel. Typically used for copy number data type. Default is FALSE.}
  \item{cap}{Image color option}
}
\seealso{
  \code{\link{iClusterBayes}},\code{\link{plotHeatmap}}
}
\details{
  The samples are ordered by the cluster assignment by the R code: order(fit$clusters).
  For each data set, the features are ordered by hierarchical clustering of the features using the complete method and 1-correlation coefficient as the distance.
  }
\value{
  no value returned.
}
\examples{
# see iManual.pdf
}
\author{Ronglai Shen \email{shenr@mskcc.org},Qianxing Mo \email{qianxing.mo@moffitt.org}}
\references{
Mo Q, Shen R, Guo C, Vannucci M, Chan KS, Hilsenbeck SG. (2018).  
A fully Bayesian latent variable model for integrative clustering analysis of multi-type omics data.
Biostatistics 19(1):71-86. 
}
\keyword{models}

