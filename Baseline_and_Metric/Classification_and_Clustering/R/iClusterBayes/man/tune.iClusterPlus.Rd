\name{tune.iClusterPlus}
\alias{tune.iClusterPlus}
\title{Integrative clustering of multiple genomic data}
\description{
Given multiple genomic data (e.g., copy number, gene
expression, DNA methylation) measured in the same set of samples,
tune.iClusterPlus uses a series of lambda values to fit a
regularized latent variable model based clustering that generates an integrated
cluster assignment based on joint inference across data.
}
\usage{
tune.iClusterPlus(cpus=8,dt1,dt2=NULL,dt3=NULL,dt4=NULL,
  type=c("gaussian","binomial","poisson","multinomial"),
  K=2,alpha=c(1,1,1,1),n.lambda=NULL,scale.lambda=c(1,1,1,1),
  n.burnin=200,n.draw=200,maxiter=20,sdev=0.05,eps=1.0e-4)
}
\arguments{
  \item{cpus}{Number of CPU used for parallel computation.}
  \item{dt1}{A data matrix. The rows represent samples, and the columns
    represent genomic features.}
  \item{dt2}{A data matrix. The rows represent samples, and the
    columns represent genomic features.}
  \item{dt3}{A data matrix. The rows represent samples, and the
    columns represent genomic features.}
  \item{dt4}{A data matrix. The rows represent samples, and the
    columns represent genomic features.}
  \item{type}{data type, which can be "gaussian","binomial","poisson", and"multinomial".}
  \item{K}{The number of eigen features. Given K, the number of cluster is K+1.}
  \item{alpha}{Vector of elasticnet penalty terms.  At this version of iClusterPlus, elasticnet is not used.  Therefore, all the elements of alpha are set to 1.}
  \item{n.lambda}{Number of lambda are tuned.}
  \item{scale.lambda}{A value between (0,1); the actual lambda values
    will be scale.lambda multiplying the lambda values of the uniform design.}
  \item{n.burnin}{Number of MCMC burnin.}
  \item{n.draw}{Number of MCMC draw.}
  \item{maxiter}{Maximum iteration for the EM algorithm.}
  \item{sdev}{standard deviation of random walk proposal.}
  \item{eps}{EM algorithm convergence criterion.}
}
\seealso{
 \code{\link{plotiCluster}},\code{\link{iClusterPlus}},\code{\link{iCluster2}},\code{\link{iCluster}}, \code{\link{compute.pod}}
}
\value{
  A list with the two elements 'fit' and 'lambda', where fit itself is
  a list and lambda is a matrix.  Each row of lambda is the lambda
  values used to fit iClusterPlus model.  Each component of fit is an object
  return by iClusterPlus, one-to-one corresponding to the row of lambda.
  Each component of fit has the following objects.   
  \item{alpha}{Intercept parameter for the genomic features.}
  \item{beta}{Information parameter for the genomic features.  The rows and the columns represent the genomic features and the coefficients for the latent variable, respectively.}
  \item{clusters}{Cluster assignment.}
  \item{centers}{Cluster centers.}
  \item{meanZ}{Latent variable.}
}
\examples{
### see the users' guide iManul.pdf 
}
\author{Qianxing Mo \email{qianxing.mo@moffitt.org}, Ronglai Shen \email{shenr@mskcc.org}}
\references{
	Qianxing Mo, Sijian Wang, Venkatraman E. Seshan, Adam B. Olshen, Nikolaus Schultz, Chris Sander, R. Scott Powers, Marc Ladanyi, and Ronglai Shen.  (2012).  
	Pattern discovery and cancer gene identification in integrated cancer genomic data. 
	Proc. Natl. Acad. Sci. USA 110(11):4245-50. 
}
\keyword{models}

