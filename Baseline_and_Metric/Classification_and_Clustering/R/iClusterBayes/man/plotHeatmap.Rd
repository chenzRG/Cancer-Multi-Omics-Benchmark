\name{plotHeatmap}
\alias{plotHeatmap}
\title{A function to generate heatmap panels sorted by integrated cluster assignment.}
\description{
  A function to generate heatmap panels sorted by integrated cluster assignment.
}
\usage{
plotHeatmap(fit,datasets,type=c("gaussian","binomial","poisson","multinomial"),
	sample.order=NULL,row.order=NULL,sparse=NULL,threshold=rep(0.25,length(datasets)),
	width=5,scale=rep("none",length(datasets)),col.scheme=rep(list(bluered(256)),
	length(datasets)), chr=NULL, plot.chr=NULL, cap=NULL)
}
\arguments{
  \item{fit}{A iCluster object.}
  \item{datasets}{A list object of data matrices.}
  \item{type}{Types of data in the datasets.}
  \item{sample.order}{User supplied cluster assignment.}
  \item{row.order}{A vector of logical values each specificy whether the genomic features in the corresponding data matrix should be reordered by similarity. Default is TRUE.}
  \item{sparse}{A vector of logical values each specificy whether to plot the top cluster-discriminant features. Default is FALSE.}
  \item{threshold}{When sparse is TRUE, a vector of threshold values to include the genomic features for which the absolute value of the associated coefficient estimates fall in the top quantile. threshold=c(0.25,0.25) takes the top quartile most discriminant features in data type 1 and data type 2 for plot.}
  \item{width}{Width of the figure in inches}
  \item{scale}{A vector of logical values each specify whether data should be scaled. Default is FALSE.}
  \item{col.scheme}{Color scheme. Can use bluered(n) in gplots R package.}
  \item{chr}{A vector of chromosome number.}
  \item{plot.chr}{A vector of logical values each specificy whether to annotate chromosome number on the left of the panel. Typically used for copy number data type. Default is FALSE.}
  \item{cap}{Image color option}
}
\details{
  The samples are ordered by the cluster assignment using the R code: order(fit$clusters).
  For each data set, the features are ordered by hierarchical clustering of the features using the complete method and 1-correlation coeffient as the distance.
  }
\seealso{
  \code{\link{iCluster}},\code{\link{iCluster2}}
}
\value{
  no value returned.
}
\examples{
# see iManual.pdf
}
\author{Ronglai Shen \email{shenr@mskcc.org}}
\references{
  Ronglai Shen, Adam Olshen, Marc Ladanyi. (2009). Integrative
  clustering of multiple genomic data types using a joint latent
  variable model with application to breast and lung cancer subtype
  analysis. \emph{Bioinformatics} 25, 2906-2912.
  
  Ronglai Shen, Qianxing Mo, Nikolaus Schultz, Venkatraman E. Seshan, Adam B. Olshen, Jason Huse, Marc Ladanyi, Chris Sander. (2012). Integrative Subtype Discovery in Glioblastoma Using iCluster. \emph{PLoS ONE} 7, e35236 
    
}
\keyword{models}

