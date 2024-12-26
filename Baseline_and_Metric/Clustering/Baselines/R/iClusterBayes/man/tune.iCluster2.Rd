\name{tune.iCluster2}
\alias{tune.iCluster2}
\title{Integrative clustering of multiple genomic data types}
\description{
Given multiple genomic data types (e.g., copy number, gene
expression, DNA methylation) measured in the same set of samples, iCluster fits a
regularized latent variable model based clustering that generates an integrated
cluster assignment based on joint inference across data types
}
\usage{
tune.iCluster2(x, K, method=c("lasso","enet","flasso","glasso","gflasso"),base=200,
  chr=NULL,true.class=NULL,lambda=NULL,n.lambda=NULL,save.nonsparse=F,nrep=10,eps=1e-4)
}
\arguments{
  \item{x}{A list object containing m data matrices representing m
    different genomic data types measured in a set of n samples. For each matrix, the
    rows represent samples, and the columns represent genomic features.}
  \item{K}{Number of subtypes.}
  \item{lambda}{User supplied matrix of lambda to tune.}
  \item{method}{Method used for clustering and variable selection.}
  \item{chr}{Chromosome labels}
  \item{n.lambda}{Number of lambda to sample using uniform design.}
  \item{nrep}{Fold of cross-validation.}
  \item{base}{Base.}
  \item{true.class}{True class label if available.}
  \item{save.nonsparse}{Logic argument whether to save the nonsparse fit.}
  \item{eps}{EM algorithm convergence criterion}  
}
\seealso{
 \code{\link{iCluster2}}
}
\value{
  A list with the following elements.
  \item{best.fit}{Best fit.}
  \item{best.lambda}{Best lambda.}
  \item{ps}{Rand index}
  \item{ps.adjusted}{Adjusted Rand index.}
}
\author{Qianxing Mo \email{qianxing.mo@moffitt.org},Ronglai Shen,Sijian Wang}
\references{
  Ronglai Shen, Sijian Wang, Qianxing Mo. (2013).  Sparse Integrative Clustering of Multiple Omics Data Sets. \emph{Annals of Applied Statistics.} 7(1):269-294
}
\keyword{models}

