\name{gbm}
\docType{data}
\alias{gbm.seg}
\alias{gbm.exp}
\alias{gbm.mut}
\title{GBM data}
\description{
 This is a subset of the glioblastoma dataset from the cancer genome atlas (TCGA) GBM study (2009) used in Shen et al. (2012).
}
\usage{data(gbm)}
\format{
  A list object containing three data matrices: copy number, methylation and mRNA expression in 84 samples.
} 
\value{
  \item{gbm.seg}{GBM copy number segmentation results genereated by DNAcopy package.} 
  \item{gbm.exp}{GBM gene expression data.}
  \item{gbm.mut}{GBM mutation data.}
}

\references{  
Ronglai Shen, Qianxing Mo, Nikolaus Schultz, Venkatraman E. Seshan, Adam B. Olshen, Jason Huse, Marc Ladanyi, Chris Sander. (2012). Integrative Subtype Discovery in Glioblastoma Using iCluster. \emph{PLoS ONE} 7, e35236 
}
\keyword{datasets}
