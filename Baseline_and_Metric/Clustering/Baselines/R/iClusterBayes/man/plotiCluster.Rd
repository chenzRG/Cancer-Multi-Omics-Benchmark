\name{plotiCluster}
\alias{plotiCluster}
\title{A function to generate cluster separability matrix plot.}
\description{
  A function to generate cluster separability matrix plot.
}
\usage{
plotiCluster(fit,label=NULL)
}
\arguments{
  \item{fit}{A iCluster object}
  \item{label}{Sample labels}
}
\seealso{
  \code{\link{iCluster}}, \code{\link{compute.pod}}
}
\value{
  no value returned.
}
\examples{
# library(iCluster)
# data(breast.chr17)
# fit=iCluster(datasets=breast.chr17, k=4, lambda=c(0.2,0.2))
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

