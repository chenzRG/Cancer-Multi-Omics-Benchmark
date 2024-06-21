\name{iCluster2}
\alias{iCluster2}
\title{Integrative clustering of multiple genomic data types}
\description{
Given multiple genomic data types (e.g., copy number, gene
expression, DNA methylation) measured in the same set of samples, iCluster fits a
regularized latent variable model based clustering that generates an integrated
cluster assigment based on joint inference across data types
}
\usage{
iCluster2(x, K, lambda, method=c("lasso","enet","flasso","glasso","gflasso"),
  chr=NULL, maxiter=50, eps=1e-4, eps2=1e-8)
}
\arguments{
  \item{x}{A list object containing m data matrices representing m
    different genomic data types measured in a set of n samples. For each matrix, the
    rows represent samples, and the columns represent genomic features.}
  \item{K}{Number of subtypes.}
  \item{lambda}{A list with m elements; each element is a vector with
    one or two elements depending on the methods used.}
  \item{method}{Method used for clustering and variable selection.}
  \item{chr}{Chromosome labels}
  \item{maxiter}{Maximum iteration for the EM algorithm.}
  \item{eps}{EM algorithm convegence criterion 1.}
  \item{eps2}{EM algorithm convegence criterion 2.}  
}
\seealso{
 \code{\link{plotiCluster}}, \code{\link{compute.pod}}, \code{\link{iClusterPlus}}
}
\value{
  A list with the following elements.
  \item{cluster}{Cluster assigment.}
  \item{centers}{cluster centers.}
  \item{Phivec}{parameter phi; a vector.}
  \item{beta}{parameter B; a matrix.}
  \item{meanZ}{meanZ}
  \item{EZZt}{EZZt}
  \item{dif}{difference}
  \item{iter}{iteration}
}
\examples{

## clustering
n1 = 20
n2 = 20
n3 = 20
n = n1+n2+n3
p = 5
q = 100

x = NULL
x1a = matrix(rnorm(n1*p), ncol=p)
x2a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
x3a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[1]] = cbind(xa,xb)

x1a = matrix(rnorm(n1*p), ncol=p)
x2a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
x3a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[2]] = cbind(xa,xb)

x1a = matrix(rnorm(n1*p), ncol=p)
x2a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
x3a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[3]] = cbind(xa,xb)


x1a = matrix(rnorm(n1*p), ncol=p)
x2a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
x3a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[4]] = cbind(xa,xb)

x1a = matrix(rnorm(n1*p), ncol=p)
x2a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
x3a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[5]] = cbind(xa,xb)

method = c('lasso', 'enet', 'flasso', 'glasso', 'gflasso')  
lambda=alist()
lambda[[1]] = 30
lambda[[2]] = c(20,1)
lambda[[3]] = c(20,20)
lambda[[4]] = 30
lambda[[5]] = c(30,20)

chr=c(rep(1,10),rep(2,(p+q)-10)) 
date()
fit2 = iCluster2(x, K=3, lambda, method=method, chr=chr, maxiter=20,eps=1e-4, eps2=1e-8)
date()

par(mfrow=c(5,1),mar=c(4,4,1,1))
for(i in 1:5){
barplot(fit2$beta[[i]][,1])
}

#library(gplots)
#library(lattice)

#plotHeatmap(fit2, datasets=x, type=rep("gaussian",length(x)),
   #row.order=c(TRUE,TRUE,FALSE,TRUE,FALSE), 
   #sparse=rep(FALSE,length(x)), scale=rep("row",5), width=5,
   #col.scheme=rep(list(bluered(256)),length(x)))

}
\author{Qianxing Mo \email{qianxing.mo@moffitt.org},Ronglai Shen,Sijian Wang}
\references{
  Ronglai Shen, Sijian Wang, Qianxing Mo. (2013).  Sparse Integrative Clustering of Multiple Omics Data Sets. \emph{Annals of Applied Statistics.} 7(1):269-294
}
\keyword{models}

