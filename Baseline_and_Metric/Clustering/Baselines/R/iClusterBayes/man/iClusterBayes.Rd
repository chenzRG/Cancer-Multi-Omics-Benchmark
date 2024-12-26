\name{iClusterBayes}
\alias{iClusterBayes}
\title{Integrative clustering of multiple genomic data types}
\description{
Given multiple genomic data types (e.g., copy number, gene
expression, DNA methylation) measured in the same set of samples, iClusterBayes fits a
Bayesian latent variable model that generates an integrated
cluster assignment based on joint inference across data types and identifies genomic features
that contribute to the clusters. 
}
\usage{
iClusterBayes(dt1,dt2=NULL,dt3=NULL,dt4=NULL,dt5=NULL,dt6=NULL,
	type = c("gaussian","binomial","poisson"),K=2,n.burnin=1000,n.draw=1200,
	prior.gamma=rep(0.1,6),sdev=0.5,beta.var.scale=1,thin=1,pp.cutoff=0.5)
}
\arguments{
  \item{dt1}{Data set 1 - a matrix with rows and columns representing samples and genomic features, respectively.}
  \item{dt2}{Data set 2 - a matrix with rows and columns representing samples and genomic features, respectively.}
  \item{dt3}{Data set 3 - a matrix with rows and columns representing samples and genomic features, respectively.}
  \item{dt4}{Data set 4 - a matrix with rows and columns representing samples and genomic features, respectively.}	    
  \item{dt5}{Data set 5 - a matrix with rows and columns representing samples and genomic features, respectively.}
  \item{dt6}{Data set 6 - a matrix with rows and columns representing samples and genomic features, respectively.}	  
  \item{type}{Data type corresponding to dt1-6, which can be gaussian, binomial, or poisson.}
  \item{K}{The number of eigen features. Given K, the number of cluster is K+1.}
  \item{n.burnin}{Number of MCMC burnin.}
  \item{n.draw}{Number of MCMC draw.}
  \item{prior.gamma}{Prior probability for the indicator variable gamma of each data set.}
  \item{sdev}{Standard deviation of random walk proposal for the latent variable.}
  \item{beta.var.scale}{A positive value to control the scale of covariance matrix of the proposed beta.}
  \item{thin}{A parameter to thin the MCMC chain in order to reduce autocorrelation.
  		Discard all but every 'thin'th sampling values. When thin=1, all sampling values are kept.}
  \item{pp.cutoff}{Posterior probability cutoff for the indicator variable gamma.  The BIC and deviance ratio will be
  	calculated by setting parameter beta to zero when the posterior probability of gamma <= cutoff.}
}
\seealso{
 \code{\link{tune.iClusterBayes}},\code{\link{plotHMBayes}},\code{\link{iClusterPlus}},\code{\link{tune.iClusterPlus}},\code{\link{plotHeatmap}}
}
\value{
  A list with the following elements.
  \item{alpha}{Intercept parameter.}
  \item{beta}{Information parameter.}
  \item{beta.pp}{Posterior probability of beta.  The higher the beta.pp, the more likely the beta should be included in the model.}
  \item{gamma.ar}{Acceptance ratio for the parameter gamma.}
  \item{beta.ar}{Acceptance ratio for the parameter beta.}
  \item{Z.ar}{Acceptance ratio for the latent variable.}
  \item{clusters}{Cluster assignment.}
  \item{centers}{Cluster center.}
  \item{meanZ}{The latent variable.}
  \item{BIC}{Bayesian information criterion.}
  \item{dev.ratio}{ see dev.ratio defined in glmnet package.}
}
\examples{
# see iManual.pdf
}
\author{Qianxing Mo \email{qianxing.mo@moffitt.org}}
\references{
Mo Q, Shen R, Guo C, Vannucci M, Chan KS, Hilsenbeck SG. (2018).  
A fully Bayesian latent variable model for integrative clustering analysis of multi-type omics data.
Biostatistics 19(1):71-86. 
}
\keyword{models}
