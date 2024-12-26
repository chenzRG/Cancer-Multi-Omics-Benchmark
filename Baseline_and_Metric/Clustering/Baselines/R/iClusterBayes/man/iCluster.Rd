\name{iCluster}
\alias{iCluster}
\title{Integrative clustering of multiple genomic data types}
\description{
Given multiple genomic data types (e.g., copy number, gene
expression, DNA methylation) measured in the same set of samples, iCluster fits a
regularized latent variable model based clustering that generates an integrated
cluster assigment based on joint inference across data types
}
\usage{
iCluster(datasets, k, lambda, scalar=FALSE, max.iter=50,epsilon=1e-3)
}
\arguments{
  \item{datasets}{A list object containing m data matrices representing m
    different genomic data types measured in a set of n samples. For each matrix, the
    rows represent samples, and the columns represent genomic features.}
  \item{k}{Number of subtypes.}
  \item{lambda}{Vector of length-m lasso penalty terms.}
  \item{scalar}{If TRUE, assumes scalar covariance matrix Psi. Default
    is FALSE.}
  \item{max.iter}{Maximum iteration for the EM algorithm.}
  \item{epsilon}{EM algorithm convegence criterion.}
}
\seealso{
 \code{\link{breast.chr17}},\code{\link{plotiCluster}}, \code{\link{compute.pod}}
}
\value{
  A list with the following elements.
  \item{meanZ}{Relaxed cluster indicator matrix.}
  \item{beta}{Coefficient matrix.}
  \item{clusters}{Cluster assigment.}
  \item{conv.rate}{Convergence history.}
}
\examples{

data(breast.chr17)
fit=iCluster(breast.chr17, k=4, lambda=c(0.2,0.2))
plotiCluster(fit=fit, label=rownames(breast.chr17[[2]]))
compute.pod(fit)

#library(gplots)
#library(lattice)
#col.scheme = alist()
#col.scheme[[1]] = bluered(256)
#col.scheme[[2]] = greenred(256)
#cn.image=breast.chr17[[2]]
#cn.image[cn.image>1.5]=1.5
#cn.image[cn.image< -1.5]= -1.5
#exp.image=breast.chr17[[1]]
#exp.image[exp.image>3]=3
#exp.image[exp.image< -3]=3
#plotHeatmap(fit, datasets=list(cn.image,exp.image), type=c("gaussian","gaussian"),
#  row.order=c(FALSE,FALSE), width=5, col.scheme=col.scheme)
}
\author{Ronglai Shen \email{shenr@mskcc.org}}
\references{
  Ronglai Shen, Adam Olshen, Marc Ladanyi. (2009). Integrative
  clustering of multiple genomic data types using a joint latent
  variable model with application to breast and lung cancer subtype
  analysis. \emph{Bioinformatics} 25, 2906-2912.
}
\keyword{models}

