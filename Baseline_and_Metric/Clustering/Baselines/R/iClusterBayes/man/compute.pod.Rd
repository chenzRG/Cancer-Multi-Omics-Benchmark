\name{compute.pod}
\alias{compute.pod}
\title{A function to compute the proportion of deviation from perfect block diagonal matrix}
\description{
  A function to compute the proportion of deviation from perfect block diagonal matrix.
}
\usage{
compute.pod(fit)
}
\arguments{
  \item{fit}{A iCluster object}
}
\seealso{
  \code{\link{iCluster}},\code{\link{iCluster2}}, \code{\link{plotiCluster}}
}
\value{
  \item{pod}{proportion of deviation from perfect block diagonal matrix}
}
\examples{
# library(iCluster)
# data(breast.chr17)
# fit=iCluster(breast.chr17, k=4, lambda=c(0.2,0.2))
# plotiCluster(fit=fit, label=rownames(breast.chr17[[2]]))
# compute.pod(fit)
}
\author{Ronglai Shen \email{shenr@mskcc.org}}
\references{
  Ronglai Shen, Adam Olshen, Marc Ladanyi. (2009). Integrative
  clustering of multiple genomic data types using a joint latent
  variable model with application to breast and lung cancer subtype
  analysis. \emph{Bioinformatics} 25, 2906-2912.
}
\keyword{models}

