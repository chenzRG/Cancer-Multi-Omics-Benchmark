###############################################################
#iCluster
#last update: 03/30/2010
#
#ARGUMENTS:
#datasets: multiple genomic platform (MGP) data in list format.
#Arrange each data set such that rows are samples and columns are genes.
#k: number of clusters.
#n: number of samples (must be the same for all data types
#lambda: a length-m vector of lasso penalty parameters
#max.iter: maximum number of iteration in the EM algorithm.
#epsilon: convergence criterion.
#scalar: a logical argument specifying whether scalar covariance matrix should be used
#VALUE:
#expZ: MLE estimate of E[Z*|X] of dimension K-1 by n where n is sample size.
#clusters: cluster assignments.
################################################################


iCluster = function(datasets, k, lambda, scalar=FALSE, max.iter=50, epsilon=1e-3){

  n=dim(datasets[[1]])[1]
  m=length(datasets)
  p=unlist(lapply(1:m,function(l){dim(datasets[[l]])[2]}))
  sum.p=sum(p)
  #error checks
  if(is.list(datasets)==F)stop("datasets must be a list")
  if(sum(is.na(datasets))>0)stop("data cannot have NAs. please exclude or impute missing values.")
  if(missing(k))stop("must specify the number of sample clusters k")
  if(missing(lambda))stop("must supply lambda values")


  #create a big matrix combining multiple data types
  stacked=matrix(NA, nrow=n, ncol=sum.p)
  if (m > 1) {
    for (i in 1:m) {
      if(i==1){idx=1:p[i]}else{idx=(sum(p[1:(i-1)])+1):sum(p[1:i])}
      stacked[,idx] <- datasets[[i]]
    }
  }
  stacked=scale(stacked,center=T,scale=F)

  #initialize W_sum.px(k-1) and PSI_sum.pxsum.p using PCA solution
  svd.x=svd(stacked)
  v=svd.x$v[,1:(k-1),drop=F]
  d=svd.x$d^2
  sigma2=sum(svd.x$d[k:n])/(n-k+1)
  cc=sqrt(d[1:(k-1)]-sigma2)
  n.cc=length(cc)
  W=v%*%diag(cc,nrow=n.cc,ncol=n.cc)

  PSI=rep(sigma2,sum.p)

  sigma=W%*%t(W)
  diag.idx=1+0:(sum.p-1)*(sum.p+1)
  sigma[diag.idx]=sigma[diag.idx]+PSI
  inv.sigma=solve(sigma)

  expZ=matrix(0,nrow=k-1,ncol=n)

  iter=1
  conv.rate=1
  save.conv.rate=NULL
  cat(paste("K=",k,":",sep=""))
  while((iter<max.iter)&(conv.rate>epsilon)){
  cat(iter)

  #E-step
  expZ.old=expZ
  expZ=t(W)%*%inv.sigma%*%t(stacked)   #kxn
  varZ=diag(c(k-1),nrow=k-1,ncol=k-1)-t(W)%*%inv.sigma%*%W
  expZZ=varZ+expZ%*%t(expZ)
  inv.expZZ=solve(expZZ)

  #M-step
  W=t(stacked)%*%t(expZ)%*%inv.expZZ
  for(i in 1:m){
      if(i==1){idx=1:p[i]}else{idx=(sum(p[1:(i-1)])+1):sum(p[1:i])}
      W[idx, ]=soft(W[idx,,drop=F],lambda[i])
  }
  normW=sqrt(apply(W,2,function(x){sum(x^2)}))
  normW[normW==0]=1
  W=t(t(W)/normW)

  a=diag(t(stacked)%*%stacked-W%*%expZ%*%stacked)
  if(scalar){PSI=diag(rep(sum(a)/n/sum.p,length(a)))}else{
      PSI=a/n
  }

  sigma=W%*%t(W)
  sigma[diag.idx]=sigma[diag.idx]+PSI
  inv.sigma=solve(sigma)

  #monitor convergence based on expZ
  conv.rate=max(abs(expZ-expZ.old))
  save.conv.rate[iter]=conv.rate

  iter=iter+1

  }
  clusters=kmeans(t(expZ),k,nstart=100)$cluster
  W.list=split(data.frame(W),f=rep(1:m,each=nrow(W)/m))
  output=list(meanZ=expZ, beta=W.list, PSI=PSI, clusters=clusters, conv.rate=save.conv.rate)
  return(output)
}

###soft thresholding###
soft=function (a, para){
  b <- sort(abs(a))
  b <- abs(a) - para
  b <- (b + abs(b))/2 #|b|+ take the positive value
  b <- sign(a) * b
  b
}
