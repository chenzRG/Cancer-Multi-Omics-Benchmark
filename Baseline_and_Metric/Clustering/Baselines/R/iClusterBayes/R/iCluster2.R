# iCluster function based on the iCluster C programs
# Last Updated: 10/11/2010 

iCluster2 = function(x, K, lambda, method=c("lasso","enet","flasso","glasso","gflasso"),
  chr=NULL, maxiter=50, eps=1e-4, eps2=1e-8){
  ## x (small) is the array, x[[k]]=...
  ## X (big) is the pool-together matrix
  ## endpoints is the endpoints of each chromoson chains (for fused lasso)
  ## for example, if endpoints=c(10,20,30), it means there are three chains: x1:x10, x11:x20, x21:x30
  ## if there is only one chain, endpoints=p
  ## dim(X)=n*p
  ## X = (X1,...,Xn)^T, Xi is a n*1 column vector of the observation for subject i
  ## Z = (Z1,...,Zn)^T, Zi is a K*1 column vector of the latent variable for subject i
  ## B = (B1,...,Bp)^T, Bj is a K*1 column vector of coefficient for variable j
  
  ## standardize data
  T = length(x)
  endpoints = 0
  ID = 1
  lenID = 1
  lbda = rep(10,8)
  allMethod = c("lasso","enet","flasso","glasso","gflasso")
  methodN = rep(NA,T)
  pvec = rep(0,T)
  for(t in 1:T){
    pvec[t] = dim(x[[t]])[2]
    methodN[t] = pmatch(method[t],allMethod)
    if(is.na(methodN[t])){
      stop("Methods are not correctly specified!\n")
    }

    if((methodN[t] == 3) || (methodN[t] == 5)){
      if(is.null(chr)){
        stop("User must supply chromosome indicator when using flasso or gflasso")
      }
      endpoints=cumsum(table(chr))
      ID = setdiff((1:pvec[t]), endpoints)
      lenID = length(ID)
    }

    if(methodN[t] == 1){
      if(length(lambda[[t]])==1){
        lbda[1] = lambda[[t]]
      }else{
        stop("Error in lambda: lambda for lasso must be single value !\n")
      }
    }else if(methodN[t] == 2){
      if(length(lambda[[t]])==2){
        lbda[2:3] = lambda[[t]]
      }else{
        stop("Error in lambda: lambda for enet must be a vector with two elements!\n")
      }
    }else if(methodN[t] == 3){
      if(length(lambda[[t]])==2){ 
        lbda[4:5] = lambda[[t]]
      }else{
        stop("Error in lambda: lambda for flasso must be a vector with two elements!\n")
      }
    }else if(methodN[t] == 4){
      if(length(lambda[[t]])==1){
        lbda[6] = lambda[[t]]
      }else{
        stop("Error in lambda: lambda for glasso must be a single value!\n")
      }
    }else if(methodN[t] == 5){
      if(length(lambda[[t]])==2){ 
        lbda[7:8] = lambda[[t]]
      }else{
        stop("Error in lambda: lambda for gflasso must be a vector with two elements!\n")
      }
    }
  }

  meanx = NULL
  normx = NULL
  for (t in 1:T) {
    meanx[[t]] = apply(x[[t]],2,mean)
    #normx[[t]] = apply(x[[t]],2,sd)
    x[[t]] = scale(x[[t]],meanx[[t]],scale=F)
  }
  
  X = NULL
  for (t in 1:T) {
    X = cbind(X, x[[t]])		
  }
  n = dim(X)[1]
  p = dim(X)[2]
  K = K-1 ## due to the contrast
  XtX = t(X)%*%X
  XtXdiag = diag(XtX)
  
  ## Initial value using SVD
  svdfit = svd(X)
  B = (as.matrix(svdfit$v[,1:K])) ## use SVD for B
  Phivec = rep(1,p)
  
  ## EM
  iter = 0
  dif = 1
  
  EZ = matrix(0,nrow=K,ncol=n)
  EZZt = matrix(0,K,K)

#  dyn.unload("iCluster.so")
#  if(!is.loaded("iCluster.so")){
#    dyn.load("iCluster.so")
#  }

  ires = .C("iClusterCore",as.integer(p),as.integer(K),as.integer(n), as.double(XtXdiag),as.double(X),
    B=as.double(B),EZ=as.double(EZ),EZZt=as.double(EZZt),Phivec=as.double(Phivec),dif=as.double(dif),
    iter = as.integer(iter),as.integer(pvec),as.double(lbda),as.double(eps),as.double(eps2),
    as.integer(maxiter),as.integer(T),as.integer(methodN),as.integer(ID),as.integer(lenID),PACKAGE="iClusterPlus")
        
  Bmat = matrix(ires$B,ncol=K)
  #Bmat[Bmat==eps2]=0
  Bmat[abs(Bmat)<= 1/n]=0

#  Blist = NULL
#  s = 0
#  for (t in 1:T) {
#    Blist[[t]] = Bmat[(s+1):(s+pvec[t]),]
#    s = s + pvec[t]
#  }
  EZ = matrix(ires$EZ,ncol=n)
  EZZt = matrix(ires$EZZt,ncol=K)
 
 #if(K==2){clusters=ifelse(EZ>=0,2,1);centers=tapply(EZ,clusters,mean)}else{
  kmeans.fit=kmeans(t(EZ),K+1,nstart=100)
  clusters=kmeans.fit$cluster
  centers=kmeans.fit$centers
  Bmat.list=split(data.frame(Bmat),f=rep(1:T,unlist(lapply(1:T,function(m)ncol(x[[m]])))))
  
  return(list(clusters=clusters, centers=centers, Phivec=ires$Phivec, beta=Bmat.list, meanZ=EZ, EZZt=EZZt,
              dif=ires$dif, iter=ires$iter))
}


