\name{utility}
\alias{getBIC}
\alias{getDevR}
\alias{getClusters}
\alias{iManual}
\title{Utility functions for iClusterPlus package}
\description{
Some utility functions for processing the results produced by iClusterPlus methods.
}
\usage{
getBIC(resultList)
getDevR(resultList)
getClusters(resultList)
iManual(view=TRUE)
}
\arguments{
  \item{resultList}{A list object as shown in the following example.}
  \item{view}{A logical value TRUE or FALSE}
}
\seealso{
 \code{\link{tune.iClusterPlus}}, \code{\link{iClusterPlus}}, \code{\link{iCluster2}}
}
\value{
  \item{getBIC}{ produce a matrix containing the BIC value for each lambda and K;
    the rows correspond to the lambda (vector) and the columns correspond to the
    K latent variables.} 
  \item{getDevR}{produce a matrix containing the deviance ratio for each lambda and K;
    the rows correspond to the lambda (vector) and the columns correspond to the
    K latent variables.}
  \item{getClusters}{produce a matrix containing the cluster assigments
    for the samples under each K; the rows correspond to the samples;
    the columns correspond to the K latent variables.}
  \item{iManual}{Open the iClusterPlus User's Guide.}
}

\examples{
### see the users' guide iManual.pdf 

#data(simuResult)
#BIC = getBIC(simuResult)
#devR = getDevR(simuResult)
#clusters = getClusters(simuResult)
}
\author{Qianxing Mo \email{qianxing.mo@moffitt.org}}
\references{
  Qianxing Mo, Sijian Wang, Venkatraman E. Seshan, Adam B. Olshen, Nikolaus Schultz, Chris Sander, R. Scott Powers, Marc Ladanyi, and Ronglai Shen.  (2012).  
  Pattern discovery and cancer gene identification in integrated cancer genomic data. 
  Proc. Natl. Acad. Sci. USA (invited revision).
}
\keyword{models}

