\name{moGap}
\alias{moGap}

\title{
  Gap statistic for clustering latent variables in \code{\link{moa-class}}.
}
\description{
  Gap statitistic is a measurement of goodness of clustering result. This is a 
  convenient function to calculate the gap statistic of clustering "moa". 
}
\usage{
moGap(x, K.max, B = 100, cluster = c("kmeans", "hclust"), plot = TRUE, 
  dist.method = "euclidean", dist.diag = FALSE, dist.upper = FALSE, dist.p = 2, 
  hcl.method = "complete", hcl.members = NULL, 
  km.iter.max = 10, km.nstart = 10, 
  km.algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"), km.trace = FALSE)
}

\arguments{
  \item{x}{
  An object of class \code{\link{moa-class}} returned by \code{\link{mbpca}}.
}
  \item{K.max}{
  The maximum number of clusters to consider, passed to \code{clusGap}
}
  \item{B}{
  The number of bootstrap, passed to \code{clusGap}
}
  \item{cluster}{
  A charater string could be either "kmeans" or "hclust" to specify the 
  clustering algorithm.
}
  \item{plot}{
  Logical; whether return the gap statistic plot.
}
  \item{dist.method}{
  Distance meaurement, passed to function \code{"dist"}.
}
  \item{dist.diag}{
  Passed to function \code{"dist"}.
}
  \item{dist.upper}{
  Passed to function \code{"dist"}.
}
  \item{dist.p}{
  Passed to function \code{"dist"}.
}
  \item{hcl.method}{
  Hierarchical clustering method, passed to \code{"hclust"}
}
  \item{hcl.members}{
  Passed to \code{"hclust"}
}
  \item{km.iter.max}{
  Maximum number of iteration in kmeans, passed to \code{"kmeans"}.
}
  \item{km.nstart}{
  An integer to specify how many random sets should be chosen. passed to \code{"kmeans"}.
}
  \item{km.algorithm}{
  Kmeans algorithm, passed to \code{"kmeans"}.
}
  \item{km.trace}{
  See function \code{"kmeans"}.
}
}

\value{
It returns a list consists of five components:

"Tab", "n", "B", "FUNcluster" - see \code{clusGap}

"nClust" - the estimated number of clusters using different method, see \code{maxSE}

}

\references{
Tibshirani, R., Walther, G. and Hastie, T. (2001). 
Estimating the number of data clusters via the Gap statistic. 
Journal of the Royal Statistical Society B, 63, 411-423.

Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., Hornik, K.(2015).  
cluster: Cluster Analysis Basics and Extensions. 
R package version 2.0.1.
}

\author{
Chen Meng
}

\seealso{
Function "clusGap" in "cluster" package
Function "dist", "hclust", "kmeans"
}

\examples{
# see examples in \code{\link{mbpca}}


data("NCI60_4arrays")
moa <- mbpca(NCI60_4arrays, ncomp = 10, k = "all", method = "globalScore", option = "lambda1", 
             center=TRUE, scale=FALSE)
gap <- moGap(moa, K.max = 12, cluster = "hcl")

genes <- moaCoef(moa)
scr <- moaScore(moa)

moa2 <- moa(NCI60_4arrays, proc.row="center_ssq1", w.data="inertia", w.row=NULL, statis=FALSE)
gap2 <- moGap(moa, K.max = 12, cluster = "hcl")
}

\keyword{ gap statistic }
\keyword{ moa }
