\name{CNregions}
\alias{CNregions}
\title{A function to remove redundant copy number regions}
\description{
This function is used to reduce copy number regions. 
}
\usage{
CNregions(seg, epsilon=0.005, adaptive=FALSE, rmCNV=FALSE, cnv=NULL,
	       frac.overlap=0.5, rmSmallseg=TRUE, nProbes=15)
}
\arguments{
  \item{seg}{DNAcopy CBS segmentation output.}
  \item{epsilon}{the maximum Euclidean distance between adjacent probes tolerated for denying a nonredundant region.
  epsilon=0 is equivalent to taking the union of all unique break points across the n samples.}
  \item{adaptive}{Vector of length-m lasso penalty terms.}
  \item{rmCNV}{If TRUE, remove germline CNV.}
  \item{cnv}{A data frame containing germline CNV data.}
  \item{frac.overlap}{A parameter needed to be explain.}
  \item{rmSmallseg}{If TRUE, remove small segment.}
  \item{nProbes}{The segment length threshold below which the segment will be removed if rmSmallseq = TRUE.}
}
\seealso{
 \code{\link{breast.chr17}},\code{\link{plotiCluster}}, \code{\link{compute.pod}},\code{\link{iCluster}},\code{\link{iClusterPlus}}
}
\value{
A matrix with reduced copy number regions.
}
\examples{

#data(gbm)
#library(GenomicRanges)
#library(cluster)
#reducedM=CNregions(seg,epsilon=0,adaptive=FALSE,rmCNV=TRUE,cnv=NULL,
#  frac.overlap=0.5, rmSmallseg=TRUE,nProbes=5)

}
\author{Ronglai Shen \email{shenr@mskcc.org}}
\references{
Qianxing Mo, Sijian Wang, Venkatraman E. Seshan, Adam B. Olshen, Nikolaus Schultz, Chris Sander, R. Scott Powers, Marc Ladanyi, and Ronglai Shen.  (2013).  
Pattern discovery and cancer gene identification in integrated cancer genomic data. 
Proc. Natl. Acad. Sci. USA.
}
\keyword{models}

