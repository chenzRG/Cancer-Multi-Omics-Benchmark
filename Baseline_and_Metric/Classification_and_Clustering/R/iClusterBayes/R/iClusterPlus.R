#nohup R CMD BATCH --no-save --no-restore mixPCAchr17.R mixPCAchr17.Rout &
#Ctrl x Ctrl q change buffer status
# derived from mixPCAcore.R,
# last update: 8/11/2011

# changed warm start to cold start
# last update: 12/11/2012



dataType <- function(x,type,K){
  n = dim(x)[1]
  p = dim(x)[2]
  sigma2 = 0
  C = 0
  class = 0
  Beta = 0
  Alpha = 0
  ty = 0
  cat = 0
  con = 0
  nclass = 0

  numclass = function(y){
    length(unique(y))
  }

  if(type=="gaussian"){
    ty = 1
    con = x
    ## Initial value using SVD
    #Warm Start
    #svdfit = svd(x)
    #Beta = (as.matrix(svdfit$v[,1:K])) ## use SVD for B
    #Cold Start
    Beta=matrix(0,nrow=ncol(x),ncol=K)
    Alpha = apply(x,2,mean)
    #residue = x - svdfit$u[,1:K]%*%t(Beta) - Alpha
    residue = x - Alpha
    sigma2 = apply(residue,2,function(x){sum(x^2)/n})
  }else if(type == "binomial"){
    ty = 2
    cat = matrix(as.numeric(as.factor(x))-1,nrow=n,ncol=p)
    nclass = apply(cat,2,numclass)
    if(any(nclass != 2)){
      stop("Error: some columns of binomial data are made of categories not equal to 2, which must be removed.\n")
    }
                                        #
    #svdfit = svd(cat)
    #Warm Start
    #Beta = (as.matrix(svdfit$v[,1:K])) ## use SVD for B	
    #Cold Start
    Beta=matrix(0,nrow=ncol(x),ncol=K)
    meanx = apply(cat,2,mean)
    Alpha = log((0.01+meanx/(1-meanx+0.01)))
    nclass = rep(2,p)
  }else if(type == "poisson"){
    ty = 3
    cat = x
    #svdfit = svd(x)
    #Beta = (as.matrix(svdfit$v[,1:K])) ## use SVD for B	
    Beta=matrix(0,nrow=ncol(x),ncol=K)
    meanx = apply(x,2,mean)      #potential problem if meanx == 0
    Alpha = log(meanx)
  }else if(type == "multinomial"){
    ty = 4
    cat = matrix(NA,nrow=n,ncol=p)
#    Upc = matrix(0,p,C)
    for(i in 1:p){
      cat[,i] = as.numeric(as.factor(x[,i])) - 1
      if(length(unique(cat[,i]))==1){
        stop("Error: some columns of multinomial data are made of a single category, which must be removed.\n")
      }
#      tmp =  table(sort(x[,i]))
#      if(length(tmp)==C){
#        Upc[i,] = as.vector(tmp)
#      }else{
#        Upc[i,match(names(tmp),class)] =  as.vector(tmp)
#      }
#      Upc[i,1:nclass[i]] = as.vector(table(sort(x[,i])))
    }
    
    nclass = apply(cat,2,numclass)
    class = sort(unique(as.vector(cat)))
    C = length(class)    
    #svdfit = svd(cat)
    #beta = (as.matrix(svdfit$v[,1:K])) ## use SVD for beta pxk	
    beta=matrix(0,nrow=ncol(x),ncol=K)
    
    Alpha = matrix(NA,nrow=p,ncol=C)
    Beta = beta
    for(i in 1:C){
      Alpha[,i] = 1.0/C
      if(i>1){
        Beta = cbind(Beta,beta)
      }
    }
  }else{
    stop("Error: type is undefined! type is one of c('gaussian','binomial','poisson','multinomial').")
  }
   #list(n=n,p=p,C=C,class=class,nclass=nclass,Alpha=Alpha,Beta=Beta,sigma2=sigma2,type=ty,con=con,cat=cat)
   list(n=n,p=p,C=C,class=class,nclass=nclass,Alpha=Alpha,Beta=Beta,sigma2=sigma2,gamma=rep(1,p),Ratio=rep(0,p),acsGamma=rep(0,p),acsBeta=rep(0,p),type=ty,con=con,cat=cat)
}

#U = p x c
#beta = p  x c x k x
#alpha = p x c
#zi = 1 x k

mcmcMix <- function(data1,data2=NULL,data3=NULL,data4=NULL,ndt,sdev=0.05,initZ,n,K,n.burnin,n.draw){

  if(missing(data1)){
    stop("Error: data1 is missing \n")
  }

  ty = rep(1,4)
  p = rep(1,4)
  C = rep(1,4)
  a = as.list(1:4)
  b = as.list(1:4)
  con = as.list(1:4)
  cat = as.list(1:4)
  class = as.list(1:4)
  sigma2 = as.list(1:4)
  nclass = as.list(1:4)

  if(ndt>0){
    ty[1] = data1$type
    p[1] = data1$p
    C[1] = data1$C
    a[[1]] = data1$Alpha
    b[[1]] = data1$Beta
    if(data1$type == 4){
      b[[1]] = t(data1$Beta)  #Beta must be transposed for logMult function in giCluster.c
      class[[1]] = data1$class
      nclass[[1]] = data1$nclass
    }
    if(data1$type == 1){
      con[[1]] = data1$con
      sigma2[[1]] = data1$sigma2
    }else{
      cat[[1]] = data1$cat
    }
  }
  
  if(ndt>1){
    ty[2] = data2$type
    p[2] = data2$p
    C[2] = data2$C
    a[[2]] = data2$Alpha
    b[[2]] = data2$Beta
    if(data2$type == 4){
      b[[2]] = t(data2$Beta)   #Beta must be transposed for logMult function in giCluster.c
      class[[2]] = data2$class
      nclass[[2]] = data2$nclass      
    }
     
    if(data2$type == 1){
      con[[2]] = data2$con
      sigma2[[2]] = data2$sigma2
    }else{
      cat[[2]] = data2$cat
    }
  }
   
  if(ndt>2){
    ty[3] = data3$type
    p[3] = data3$p
    C[3] = data3$C
    a[[3]] = data3$Alpha
    b[[3]] = data3$Beta
    if(data3$type == 4){
      b[[3]] = t(data3$Beta)  #Beta must be transposed for logMult function in giCluster.c
      class[[3]] = data3$class
      nclass[[3]] = data3$nclass      
    }
    
    if(data3$type == 1){
      con[[3]] = data3$con
      sigma2[[3]] = data3$sigma2
    }else{
      cat[[3]] = data3$cat

    }
  }
   
  if(ndt>3){
    ty[4] = data4$type
    p[4] = data4$p
    C[4] = data4$C
    a[[4]] = data4$Alpha
    b[[4]] = data4$Beta
    if(data4$type == 4){
      b[[4]] = t(data4$Beta)  #Beta must be transposed for logMult function in giCluster.c
      class[[4]] = data4$class
      nclass[[4]] = data4$nclass      
    }
     
    if(data4$type == 1){
      con[[4]] = data4$con
      sigma2[[4]] = data4$sigma2
    }else{
      cat[[4]] = data4$cat
    }
  }
  
  meanZ = matrix(0,nrow=n,ncol=K)
  res = .C("mcmcMix",meanZ=as.double(meanZ),lastZ=as.double(initZ),as.integer(n),as.integer(K),
    as.integer(n.burnin),as.integer(n.draw),as.double(sdev),as.integer(ndt),
    as.integer(ty[1]),as.integer(p[1]),as.integer(C[1]),as.double(a[[1]]),as.double(b[[1]]),
    as.double(con[[1]]),as.integer(cat[[1]]),as.integer(class[[1]]),as.integer(nclass[[1]]),
    as.double(sigma2[[1]]),
    as.integer(ty[2]),as.integer(p[2]),as.integer(C[2]),as.double(a[[2]]),as.double(b[[2]]),
    as.double(con[[2]]),as.integer(cat[[2]]),as.integer(class[[2]]),as.integer(nclass[[2]]),
    as.double(sigma2[[2]]),
    as.integer(ty[3]),as.integer(p[3]),as.integer(C[3]),as.double(a[[3]]),as.double(b[[3]]),
    as.double(con[[3]]),as.integer(cat[[3]]),as.integer(class[[3]]),as.integer(nclass[[3]]),
    as.double(sigma2[[3]]),
    as.integer(ty[4]),as.integer(p[4]),as.integer(C[4]),as.double(a[[4]]),as.double(b[[4]]),
    as.double(con[[4]]),as.integer(cat[[4]]),as.integer(class[[4]]),as.integer(nclass[[4]]),
    as.double(sigma2[[4]]),accept=integer(n),PACKAGE="iClusterPlus")
  
 list(meanZ=matrix(res$meanZ,nrow=n,ncol=K),lastZ=matrix(res$lastZ,nrow=n,ncol=K),accept=res$accept)
}

########## total BIC for all data sets #########
totalBIC = function(Data,meanZ,ndt,K){
  BIC = 0
  loglike = 0
  for(i in 1:ndt){
    if(Data[[i]]$type == 1){  # normal #
      fit1 = .C("logNormAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha),
        as.double(Data[[i]]$Beta),as.double(Data[[i]]$sigma2),as.double(Data[[i]]$con),
        as.integer(Data[[i]]$n),as.integer(Data[[i]]$p),as.integer(K),PACKAGE="iClusterPlus")
      BIC = BIC - 2*fit1$loglike + sum(Data[[i]]$Beta != 0)*log(Data[[i]]$n * Data[[i]]$p)
      #BIC = BIC - 2*fit1$loglike + (Data[[i]]$p+sum(Data[[i]]$Beta != 0))*log(Data[[i]]$n)
    }else if(Data[[i]]$type == 2){ # binomial #
      fit2 = .C("logBinomAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha),
        as.double(Data[[i]]$Beta),as.integer(Data[[i]]$cat),
        as.integer(Data[[i]]$n),as.integer(Data[[i]]$p),as.integer(K),PACKAGE="iClusterPlus")
      BIC = BIC - 2*fit2$loglike + sum(Data[[i]]$Beta != 0)*log(Data[[i]]$n * Data[[i]]$p)
      #BIC = BIC - 2*fit2$loglike + (Data[[i]]$p+sum(Data[[i]]$Beta != 0))*log(Data[[i]]$n)
    }else if(Data[[i]]$type == 3){ # Poisson #
      fit3 = .C("logPoissonAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha),
        as.double(Data[[i]]$Beta),as.integer(Data[[i]]$cat),
        as.integer(Data[[i]]$n),as.integer(Data[[i]]$p),as.integer(K),PACKAGE="iClusterPlus")
      BIC = BIC - 2*fit3$loglike + sum(Data[[i]]$Beta != 0)*log(Data[[i]]$n * Data[[i]]$p)
      #BIC = BIC - 2*fit3$loglike + (Data[[i]]$p+sum(Data[[i]]$Beta != 0))*log(Data[[i]]$n)
    }else { # Multinomial, Beta must be t(Beta) for logMult function in giCluster.c #
      fit4 = .C("logMultAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha),
        as.double(t(Data[[i]]$Beta)),as.integer(Data[[i]]$cat),as.integer(Data[[i]]$class),
        as.integer(Data[[i]]$nclass),as.integer(Data[[i]]$n),as.integer(Data[[i]]$p),
        as.integer(Data[[i]]$C),as.integer(K),PACKAGE="iClusterPlus")
      BIC = BIC - 2*fit4$loglike + sum(Data[[i]]$Beta != 0)*log(Data[[i]]$n * Data[[i]]$p)/Data[[i]]$C
      #BIC = BIC - 2*fit4$loglike + (Data[[i]]$p+sum(Data[[i]]$Beta != 0))*log(Data[[i]]$n)
    }
  }
  BIC
}

########### deviance ratio function ############
dev.ratio = function(Data,meanZ,alpha,lambda,ndt,K){
  sumdev0 = 0
  sumdev = 0
  for(i in 1:ndt){
    if(Data[[i]]$type == 1){  # normal #
      fit1 = .C("elnetBatchDev",a0=double(Data[[i]]$p),beta=double(Data[[i]]$p*K),sigma2 = double(Data[[i]]$p),
        as.double(meanZ),as.double(Data[[i]]$con),as.integer(Data[[i]]$n),as.integer(K),as.integer(Data[[i]]$p),
        as.double(alpha[i]),as.double(lambda[i]),sumdev0=double(1),sumdev=double(1),PACKAGE="iClusterPlus")
      sumdev0 = sumdev0 + fit1$sumdev0
      sumdev = sumdev + fit1$sumdev
    }else if(Data[[i]]$type == 2){ # binomial #
      fit2 = .C("lognetBatchDev",a0=double(Data[[i]]$p),beta=double(Data[[i]]$p*K),as.double(meanZ),
        as.integer(Data[[i]]$cat),as.integer(Data[[i]]$n),as.integer(K),as.integer(Data[[i]]$p),
        as.double(alpha[i]),as.double(lambda[i]),as.integer(Data[[i]]$nclass),as.integer(2),
        as.integer(0),sumdev0=double(1),sumdev=double(1),PACKAGE="iClusterPlus") #family=0 is binomial
      sumdev0 = sumdev0 + fit2$sumdev0
      sumdev = sumdev + fit2$sumdev
    }else if(Data[[i]]$type == 3){ # Poisson #
      fit3 = .C("fishnetBatchDev",a0=double(Data[[i]]$p),beta=double(Data[[i]]$p*K),as.double(meanZ),
        as.double(Data[[i]]$cat),as.integer(Data[[i]]$n),as.integer(K),as.integer(Data[[i]]$p),
        as.double(alpha[i]),as.double(lambda[i]),sumdev0=double(1),sumdev=double(1),PACKAGE="iClusterPlus")
      sumdev0 = sumdev0 + fit3$sumdev0
      sumdev = sumdev + fit3$sumdev
    }else { # Multinomial #
      fit4 = .C("lognetBatchDev",a0=double(Data[[i]]$p * Data[[i]]$C),beta=double(Data[[i]]$p*K*Data[[i]]$C),
        as.double(meanZ),as.integer(Data[[i]]$cat),as.integer(Data[[i]]$n),as.integer(K),
        as.integer(Data[[i]]$p),as.double(alpha[i]),as.double(lambda[i]),as.integer(Data[[i]]$nclass),
        as.integer(Data[[i]]$C),as.integer(1),sumdev0=double(1),sumdev=double(1),PACKAGE="iClusterPlus") #family=1 is multinomial
      sumdev0 = sumdev0 + fit4$sumdev0
      sumdev = sumdev + fit4$sumdev
    }
  }

  # deviance ratio, for linear reg, it is R-square; the bigger, the better
  1-sumdev/sumdev0
}


iClusterPlus = function(dt1,dt2=NULL,dt3=NULL,dt4=NULL,type=c("gaussian","binomial","poisson","multinomial"),
  K=2,alpha=c(1,1,1,1),lambda=c(0.03,0.03,0.03,0.03),n.burnin=100,n.draw=200,maxiter=20,sdev=0.05,eps=1.0e-4){

  dttype = c("gaussian","binomial","poisson","multinomial")
  if(missing(dt1)){
    stop("Error: dt1 is missing!\n")
  }

  if(!all(type %in% dttype)){
      cat("Error: ",type[!all(type %in% dttype)],"\n")
      stop("Allowed data types are gaussian, binomial, poisson and multinomial. \n")
  }

  isNULL = c(is.null(dt1),is.null(dt2),is.null(dt3),is.null(dt4))
  if(any(diff(isNULL) == -1)){
      stop("Error: dt1 to dt4 must be assigned in order.\n")
  }

  if(sum(!isNULL) > length(type)){
      stop("Error:  data type is missing for some data. \n")
  } 

  n = nrow(dt1)
  ndt = 1
  Alpha = list()
  Beta = list()
  Data = list()
  Dif = list()

  Data[[1]] = dataType(dt1,type[1],K)
  Alpha[[1]] = Data[[1]]$Alpha
  Beta[[1]] = Data[[1]]$Beta
  
  Data[[2]] = NULL
  if(!is.null(dt2)){
    if(n != nrow(dt2)){stop("Error: nrow(dt1) != nrow(dt2) \n")} 
    Data[[2]] = dataType(dt2,type[2],K)
    Alpha[[2]] = Data[[2]]$Alpha
    Beta[[2]] = Data[[2]]$Beta
    ndt = ndt + 1
  }

  Data[[3]] = NULL
  if(!is.null(dt3)){
    if(n != nrow(dt3)){stop("Error: nrow(dt1) != nrow(dt3) \n")} 
    Data[[3]] = dataType(dt3,type[3],K)
    Alpha[[3]] = Data[[3]]$Alpha
    Beta[[3]] = Data[[3]]$Beta
    ndt = ndt + 1
  }

  Data[[4]] = NULL
  if(!is.null(dt4)){
    if(n != nrow(dt4)){stop("Error: nrow(dt1) != nrow(dt4) \n")} 
    Data[[4]] = dataType(dt4,type[4],K)
    Alpha[[4]] = Data[[4]]$Alpha
    Beta[[4]] = Data[[4]]$Beta
    ndt = ndt + 1
  }
  
  dif = 1
  iter = 0
  lastZ = matrix(rnorm(n*K,0,1),ncol=K)
  newZ = NULL
  while(iter<maxiter && dif>eps) {
    iter = iter + 1
#    cat(iter,' ')
    newZ = mcmcMix(Data[[1]],Data[[2]],Data[[3]],Data[[4]],ndt,sdev,lastZ,n,K,n.burnin,n.draw)
    lastZ = newZ$lastZ
#    print(table(newZ$accept))
    
    for(i in 1:ndt){
      if(Data[[i]]$type == 1){  # normal #
        fit = .C("elnetBatch",a0=double(Data[[i]]$p),beta=double(Data[[i]]$p*K),sigma2 = double(Data[[i]]$p),
          as.double(newZ$meanZ),as.double(Data[[i]]$con),as.integer(n),as.integer(K),as.integer(Data[[i]]$p),
          as.double(alpha[i]),as.double(lambda[i]),PACKAGE="iClusterPlus")

        Data[[i]]$Alpha = fit$a0
        Data[[i]]$Beta = matrix(fit$beta,ncol=K)
        Data[[i]]$sigma2 = fit$sigma2
      }else if(Data[[i]]$type == 2){ # binomial #
        fit = .C("lognetBatch",a0=double(Data[[i]]$p),beta=double(Data[[i]]$p*K),as.double(newZ$meanZ),
          as.integer(Data[[i]]$cat),as.integer(n),as.integer(K),as.integer(Data[[i]]$p),
          as.double(alpha[i]),as.double(lambda[i]),as.integer(Data[[i]]$nclass),as.integer(2),
          as.integer(0),PACKAGE="iClusterPlus") #family=0 is binomial
        
        Data[[i]]$Alpha = fit$a0
        Data[[i]]$Beta = matrix(fit$beta,ncol=K)
      }else if(Data[[i]]$type == 3){ # Poisson #
        fit = .C("fishnetBatch",a0=double(Data[[i]]$p),beta=double(Data[[i]]$p*K),as.double(newZ$meanZ),
          as.double(Data[[i]]$cat),as.integer(n),as.integer(K),as.integer(Data[[i]]$p),
          as.double(alpha[i]),as.double(lambda[i]),PACKAGE="iClusterPlus")
 
        Data[[i]]$Alpha = fit$a0
        Data[[i]]$Beta = matrix(fit$beta,ncol=K)
      }else { # Multinomial #
        fit = .C("lognetBatch",a0=double(Data[[i]]$p * Data[[i]]$C),beta=double(Data[[i]]$p*K*Data[[i]]$C),
          as.double(newZ$meanZ),as.integer(Data[[i]]$cat),as.integer(n),as.integer(K),as.integer(Data[[i]]$p),
          as.double(alpha[i]),as.double(lambda[i]),as.integer(Data[[i]]$nclass),
          as.integer(Data[[i]]$C),as.integer(1),PACKAGE="iClusterPlus") #family=1 is multinomial
        
        Data[[i]]$Alpha = matrix(fit$a0,ncol=Data[[i]]$C)
        Data[[i]]$Beta = matrix(fit$beta,ncol=K*Data[[i]]$C)       
      }
    }
 
    allDif = NULL
    for(i in 1:ndt){
      Dif[[i]] = abs(cbind((Data[[i]]$Alpha - Alpha[[i]]),(Data[[i]]$Beta - Beta[[i]])))
#      print(summary(as.vector(Dif[[i]])))
      allDif = c(allDif,as.vector(Dif[[i]]))
      Alpha[[i]] = Data[[i]]$Alpha
      Beta[[i]] = Data[[i]]$Beta
    }
    dif = max(allDif)
  }
#  cat('\n')
  for(i in 1:ndt){
    Dif[[i]] = summary(Dif[[i]])
    colnames(Dif[[i]])=c("alpha",paste("beta",1:(ncol(Dif[[i]])-1),sep=""))
  }
  
  kmeans.fit=kmeans(newZ$meanZ,K+1,nstart=100)
  clusters=kmeans.fit$cluster
  centers=kmeans.fit$centers

  BIC = totalBIC(Data,newZ$meanZ,ndt,K)
  devRatio = dev.ratio(Data,newZ$meanZ,alpha,lambda,ndt,K)
  list(alpha=Alpha,beta=Beta,clusters=clusters,centers=centers,meanZ=newZ$meanZ,
       BIC=BIC,dev.ratio=devRatio,dif=Dif)
}

