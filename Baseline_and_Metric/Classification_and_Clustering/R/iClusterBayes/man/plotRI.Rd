\name{plotRI}
\alias{plotRI}
\title{A function to generate reproducibility index plot.}
\description{
  A function to generate reproducibility index plot.
}
\usage{
plotRI(cv.fit)
}
\arguments{
  \item{cv.fit}{A tune.iCluster2 object}
}
\seealso{
  \code{\link{iCluster}}
}
\value{
  no value returned.
}
\examples{
#data(simu.datasets)
#cv.fit=alist()
#for(k in 2:5){
#  cat(paste("K=",k,sep=""),'\n')
#  cv.fit[[k]]=tune.iCluster2(datasets=simu.datasets, k,nrep=2, n.lambda=8)
#}

##Reproducibility index (RI) plot
#plotRI(cv.fit)
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

