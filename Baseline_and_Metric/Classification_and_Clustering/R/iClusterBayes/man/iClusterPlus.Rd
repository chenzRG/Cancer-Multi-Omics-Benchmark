\name{iClusterPlus}
\alias{iClusterPlus}
\title{Integrative clustering of multiple genomic data types}
\description{
Given multiple genomic data types (e.g., copy number, gene
expression, DNA methylation) measured in the same set of samples, iClusterPlus fits a
regularized latent variable model based clustering that generates an integrated
cluster assignment based on joint inference across data types
}
\usage{
iClusterPlus(dt1,dt2=NULL,dt3=NULL,dt4=NULL,
	type=c("gaussian","binomial","poisson","multinomial"),
 	K=2,alpha=c(1,1,1,1),lambda=c(0.03,0.03,0.03,0.03),
	n.burnin=100,n.draw=200,maxiter=20,sdev=0.05,eps=1.0e-4)
}
\arguments{
  \item{dt1}{A data matrix. The rows represent samples, and the columns
    represent genomic features.}
  \item{dt2}{A data matrix. The rows represent samples, and the
    columns represent genomic features.}
  \item{dt3}{A data matrix. The rows represent samples, and the
    columns represent genomic features.}
  \item{dt4}{A data matrix. The rows represent samples, and the
    columns represent genomic features.}
  \item{type}{Data type, which can be gaussian, binomial, poisson, multinomial.}
  \item{K}{The number of eigen features. Given K, the number of cluster is K+1.}
  \item{alpha}{Vector of elasticnet penalty terms. At this version of iClusterPlus, elasticnet is not used.  Therefore, all the elements of alpha are set to 1.}
  \item{lambda}{Vector of lasso penalty terms.}
  \item{n.burnin}{Number of MCMC burnin.}
  \item{n.draw}{Number of MCMC draw.}
  \item{maxiter}{Maximum iteration for the EM algorithm.}
  \item{sdev}{standard deviation of random walk proposal.}
  \item{eps}{Algorithm convergence criterion.}
}
\seealso{
 \code{\link{plotiCluster}},\code{\link{iCluster}}, \code{\link{compute.pod}}
}
\value{
  A list with the following elements.
  \item{alpha}{Intercept parameter.}
  \item{beta}{Information parameter.  }
  \item{clusters}{Cluster assignment.}
  \item{centers}{Cluster center.}
  \item{meanZ}{Latent variable.}
  \item{BIC}{Bayesian information criterion.}
  \item{dev.ratio}{ see dev.ratio defined in glmnet package.}
  \item{dif}{absolute difference for the parameters in the last and
    next-to-last iterations.}
%  \item{lastZ}{Last draw of latent variable.}
}
\examples{
# see iManual.pdf
}
\author{Qianxing Mo \email{qianxing.mo@moffitt.org},Ronglai Shen, Sijian Wang}
\references{
Qianxing Mo, Sijian Wang, Venkatraman E. Seshan, Adam B. Olshen, Nikolaus Schultz, Chris Sander, R. Scott Powers, Marc Ladanyi, and Ronglai Shen.  (2013).  
Pattern discovery and cancer gene identification in integrated cancer genomic data. 
\emph{Proc. Natl. Acad. Sci. USA.} 110(11):4245-50. 
}
\keyword{models}

