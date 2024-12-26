##### R function for iClusterBayes #########
# Qianxing Mo (qmo@bcm.edu), Dan L. Duncan Cancer Center, Baylor College of Medicine
# One Baylor Plaza, Houston, TX, 77030 
# last update 11/13/2015

########## total BIC for all data sets #########
totalBICbayes = function(Data,meanZ,ndt,K,pp.cutoff){
    BIC = 0
    for(i in 1:ndt){
    #sigID = which(Data[[i]]$Ratio >= quantile(Data[[i]]$Ratio,prob=0.5))
        sigID = which(Data[[i]]$Ratio > pp.cutoff)
        lenp = length(sigID)
        if(lenp > 0){
            if(Data[[i]]$type == 1){  # normal #
                fit1 = .C("logNormAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha[sigID]),
                    as.double(Data[[i]]$Beta[sigID,]),as.double(Data[[i]]$sigma2[sigID]),as.double(Data[[i]]$con[,sigID]),
                    as.integer(Data[[i]]$n),as.integer(lenp),as.integer(K),PACKAGE="iClusterPlus")
                BIC = BIC - 2*fit1$loglike + (lenp*(K+1))*log(Data[[i]]$n)
            }else if(Data[[i]]$type == 2){ # binomial #
                fit2 = .C("logBinomAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha[sigID]),
                    as.double(Data[[i]]$Beta[sigID,]),as.integer(Data[[i]]$cat[,sigID]),
                    as.integer(Data[[i]]$n),as.integer(lenp),as.integer(K),PACKAGE="iClusterPlus")
                BIC = BIC - 2*fit2$loglike +  (lenp*(K+1))*log(Data[[i]]$n)
            }else if(Data[[i]]$type == 3){ # Poisson #
                fit3 = .C("logPoissonAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha[sigID]),
                    as.double(Data[[i]]$Beta[sigID,]),as.integer(Data[[i]]$cat[,sigID]),
                    as.integer(Data[[i]]$n),as.integer(lenp),as.integer(K),PACKAGE="iClusterPlus")
                BIC = BIC - 2*fit3$loglike + (lenp*(K+1))*log(Data[[i]]$n)
            }
        }
        noiseID = which(Data[[i]]$Ratio <= pp.cutoff)
        lenp = length(noiseID)
        ZeroBeta = matrix(0,nrow=lenp,ncol=K)
        if(lenp > 0){
            if(Data[[i]]$type == 1){  # normal #
                fit1 = .C("logNormAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha[noiseID]),
                    as.double(ZeroBeta),as.double(Data[[i]]$sigma2[noiseID]),as.double(Data[[i]]$con[,noiseID]),
                    as.integer(Data[[i]]$n),as.integer(lenp),as.integer(K),PACKAGE="iClusterPlus")
                BIC = BIC - 2*fit1$loglike + (lenp)*log(Data[[i]]$n)
            }else if(Data[[i]]$type == 2){ # binomial #
                fit2 = .C("logBinomAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha[noiseID]),
                    as.double(ZeroBeta),as.integer(Data[[i]]$cat[,noiseID]),
                    as.integer(Data[[i]]$n),as.integer(lenp),as.integer(K),PACKAGE="iClusterPlus")
                BIC = BIC - 2*fit2$loglike +  (lenp)*log(Data[[i]]$n)
            }else if(Data[[i]]$type == 3){ # Poisson #
                fit3 = .C("logPoissonAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha[noiseID]),
                    as.double(ZeroBeta),as.integer(Data[[i]]$cat[,noiseID]),
                    as.integer(Data[[i]]$n),as.integer(lenp),as.integer(K),PACKAGE="iClusterPlus")
                BIC = BIC - 2*fit3$loglike + (lenp)*log(Data[[i]]$n)
            }
        }
    }
    BIC
}

#function to calculate deviance ratio 
dev.ratio.bayes = function(Data,meanZ,ndt,K,pp.cutoff){
    #cat("- deviance ratio -\n")
    loglike = 0
    loglikeNull = 0
    loglikeFull = 0
    for(i in 1:ndt){
        sigID = which(Data[[i]]$Ratio > pp.cutoff)
        lenp=nrow(Data[[i]]$Beta)
        ZeroBeta = matrix(0,nrow=lenp,ncol=K)
        Beta = ZeroBeta
        Beta[sigID,] = Data[[i]]$Beta[sigID,]
        #cat(lenp,"\n")
        if(Data[[i]]$type == 1){  # normal #
            fit1 = .C("logNormAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha),
                as.double(Beta),as.double(Data[[i]]$sigma2),as.double(Data[[i]]$con),
                as.integer(Data[[i]]$n),as.integer(lenp),as.integer(K),PACKAGE="iClusterPlus")
            loglike = loglike + fit1$loglike
            fit0 = .C("logNormAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha),
                as.double(ZeroBeta),as.double(Data[[i]]$sigma2),as.double(Data[[i]]$con),
                as.integer(Data[[i]]$n),as.integer(lenp),as.integer(K),PACKAGE="iClusterPlus")
            loglikeNull = loglikeNull + fit0$loglike
            fit2 = .C("logNormAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha),
                as.double(Data[[i]]$Beta),as.double(Data[[i]]$sigma2),as.double(Data[[i]]$con),
                as.integer(Data[[i]]$n),as.integer(lenp),as.integer(K),PACKAGE="iClusterPlus")
            loglikeFull = loglikeFull + fit2$loglike
        }else if(Data[[i]]$type == 2){ # binomial #
            fit1 = .C("logBinomAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha),
                as.double(Beta),as.integer(Data[[i]]$cat),
                as.integer(Data[[i]]$n),as.integer(lenp),as.integer(K),PACKAGE="iClusterPlus")
            loglike = loglike + fit1$loglike
            #cat(fit1$loglike,"\n")
            fit0 = .C("logBinomAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha),
                as.double(ZeroBeta),as.integer(Data[[i]]$cat),
                as.integer(Data[[i]]$n),as.integer(lenp),as.integer(K),PACKAGE="iClusterPlus")           
            loglikeNull = loglikeNull + fit0$loglike
            fit2 = .C("logBinomAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha),
                as.double(Data[[i]]$Beta),as.integer(Data[[i]]$cat),
                as.integer(Data[[i]]$n),as.integer(lenp),as.integer(K),PACKAGE="iClusterPlus")
            loglikeFull = loglikeFull + fit2$loglike            
            #cat(fit0$loglike,"\n")
        }else if(Data[[i]]$type == 3){ # Poisson #
            fit1 = .C("logPoissonAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha),
                as.double(Beta),as.integer(Data[[i]]$cat),
                as.integer(Data[[i]]$n),as.integer(lenp),as.integer(K),PACKAGE="iClusterPlus")
            loglike = loglike + fit1$loglike
            fit0 = .C("logPoissonAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha),
                as.double(ZeroBeta),as.integer(Data[[i]]$cat),
                as.integer(Data[[i]]$n),as.integer(lenp),as.integer(K),PACKAGE="iClusterPlus")
            loglikeNull = loglikeNull + fit0$loglike
            fit2 = .C("logPoissonAll",loglike = double(1),as.double(meanZ),as.double(Data[[i]]$Alpha),
                as.double(Data[[i]]$Beta),as.integer(Data[[i]]$cat),
                as.integer(Data[[i]]$n),as.integer(lenp),as.integer(K),PACKAGE="iClusterPlus")
            loglikeFull = loglikeFull + fit2$loglike                 
        }
    }
    #(loglike-loglikeNull)/(loglikeFull-loglikeNull)
    1-loglike/loglikeNull
}


mcmcBayes <- function(dt1,dt2=NULL,dt3=NULL,dt4=NULL,dt5=NULL,dt6=NULL,ndt,z.sdev,n,K,zBurnin,zDraw,betaBurnin,betaDraw,
                      thin,pg,beta0,invSigma0,invSigmaBeta0,invga0,invgb0,betaVarScale){
  
  if(missing(dt1)){
    stop("Error: dt1 is missing \n")
  }

  if(ndt < 1){
      stop("Error: ndt must be >= 1!")
  }

  ty = rep(0,6)
  p = rep(1,6)
  C = rep(1,6)
  a = as.list(1:6)
  b = as.list(1:6)
  con = as.list(1:6)
  cat = as.list(1:6)
  class = as.list(1:6)
  sigma2 = as.list(1:6)
  nclass = as.list(1:6)
  gamma = as.list(1:6)
  acceptBeta = as.list(1:6)
  acceptGamma = as.list(1:6)

  if(ndt>0){
    ty[1] = dt1$type
    p[1] = dt1$p
    C[1] = dt1$C
    a[[1]] = dt1$Alpha
    b[[1]] = dt1$Beta
    if(dt1$type == 4){
      b[[1]] = t(dt1$Beta)  #Beta must be transposed for logMult function in giCluster.c
      class[[1]] = dt1$class
      nclass[[1]] = dt1$nclass
    }
    if(dt1$type == 1){
        con[[1]] = dt1$con
        sigma2[[1]] = dt1$sigma2
    }else{
        cat[[1]] = dt1$cat
        sigma2[[1]] = dt1$sigma2
    }
    gamma[[1]] = dt1$gamma
    acceptBeta[[1]] = dt1$acsBeta
    acceptGamma[[1]] = dt1$acsGamma
  }
  
  if(ndt>1){
    ty[2] = dt2$type
    p[2] = dt2$p
    C[2] = dt2$C
    a[[2]] = dt2$Alpha
    b[[2]] = dt2$Beta
    if(dt2$type == 4){
      b[[2]] = t(dt2$Beta)   #Beta must be transposed for logMult function in giCluster.c
      class[[2]] = dt2$class
      nclass[[2]] = dt2$nclass      
    }
     
    if(dt2$type == 1){
      con[[2]] = dt2$con
      sigma2[[2]] = dt2$sigma2
    }else{
        cat[[2]] = dt2$cat
        sigma2[[2]] = dt2$sigma2
    }
    gamma[[2]] = dt2$gamma
    acceptBeta[[2]] = dt2$acsBeta
    acceptGamma[[2]] = dt2$acsGamma
  }
   
  if(ndt>2){
    ty[3] = dt3$type
    p[3] = dt3$p
    C[3] = dt3$C
    a[[3]] = dt3$Alpha
    b[[3]] = dt3$Beta
    if(dt3$type == 4){
      b[[3]] = t(dt3$Beta)  #Beta must be transposed for logMult function in giCluster.c
      class[[3]] = dt3$class
      nclass[[3]] = dt3$nclass      
    }
    
    if(dt3$type == 1){
      con[[3]] = dt3$con
      sigma2[[3]] = dt3$sigma2
    }else{
        cat[[3]] = dt3$cat
        sigma2[[3]] = dt3$sigma2
    }
    gamma[[3]] = dt3$gamma
    acceptBeta[[3]] = dt3$acsBeta
    acceptGamma[[3]] = dt3$acsGamma
  }
   
  if(ndt>3){
    ty[4] = dt4$type
    p[4] = dt4$p
    C[4] = dt4$C
    a[[4]] = dt4$Alpha
    b[[4]] = dt4$Beta
    if(dt4$type == 4){
      b[[4]] = t(dt4$Beta)  #Beta must be transposed for logMult function in giCluster.c
      class[[4]] = dt4$class
      nclass[[4]] = dt4$nclass      
    }
     
    if(dt4$type == 1){
      con[[4]] = dt4$con
      sigma2[[4]] = dt4$sigma2
    }else{
        cat[[4]] = dt4$cat
        sigma2[[4]] = dt4$sigma2
    }
    gamma[[4]] = dt4$gamma
    acceptBeta[[4]] = dt4$acsBeta
    acceptGamma[[4]] = dt4$acsGamma
}

  if(ndt>4){
    ty[5] = dt5$type
    p[5] = dt5$p
    C[5] = dt5$C
    a[[5]] = dt5$Alpha
    b[[5]] = dt5$Beta
    if(dt5$type == 5){
      b[[5]] = t(dt5$Beta)  #Beta must be transposed for logMult function in giCluster.c
      class[[5]] = dt5$class
      nclass[[5]] = dt5$nclass      
    }
     
    if(dt5$type == 1){
      con[[5]] = dt5$con
      sigma2[[5]] = dt5$sigma2
    }else{
        cat[[5]] = dt5$cat
        sigma2[[5]] = dt5$sigma2
    }
    gamma[[5]] = dt5$gamma
    acceptBeta[[5]] = dt5$acsBeta
    acceptGamma[[5]] = dt5$acsGamma
  }

  if(ndt>5){
    ty[6] = dt6$type
    p[6] = dt6$p
    C[6] = dt6$C
    a[[6]] = dt6$Alpha
    b[[6]] = dt6$Beta
    if(dt6$type == 6){
      b[[6]] = t(dt6$Beta)  #Beta must be transposed for logMult function in giCluster.c
      class[[6]] = dt6$class
      nclass[[6]] = dt6$nclass      
    }
     
    if(dt6$type == 1){
      con[[6]] = dt6$con
      sigma2[[6]] = dt6$sigma2
    }else{
        cat[[6]] = dt6$cat
        sigma2[[6]] = dt6$sigma2
    }
    gamma[[6]] = dt6$gamma
    acceptBeta[[6]] = dt6$acsBeta
    acceptGamma[[6]] = dt6$acsGamma
  }

  
  meanZ = matrix(0,nrow=n,ncol=K)
  initZ = matrix(rnorm(n*K,0,1),ncol=K)
  
  res = .C("mcmcBayes",sumMeanZ=as.double(meanZ),lastZ=as.double(initZ),acceptZ=as.integer(rep(0,n)),as.integer(c(n,K,zBurnin,zDraw,ndt,thin)),as.double(z.sdev),
    as.integer(c(ty[1],p[1],C[1])),alpha1=as.double(a[[1]]),beta1=as.double(b[[1]]),as.double(con[[1]]),as.integer(cat[[1]]),sigma21=as.double(sigma2[[1]]),
    as.integer(c(ty[2],p[2],C[2])),alpha2=as.double(a[[2]]),beta2=as.double(b[[2]]),as.double(con[[2]]),as.integer(cat[[2]]),sigma22=as.double(sigma2[[2]]),
    as.integer(c(ty[3],p[3],C[3])),alpha3=as.double(a[[3]]),beta3=as.double(b[[3]]),as.double(con[[3]]),as.integer(cat[[3]]),sigma23=as.double(sigma2[[3]]),
    as.integer(c(ty[4],p[4],C[4])),alpha4=as.double(a[[4]]),beta4=as.double(b[[4]]),as.double(con[[4]]),as.integer(cat[[4]]),sigma24=as.double(sigma2[[4]]),
    as.integer(c(ty[5],p[5],C[5])),alpha5=as.double(a[[5]]),beta5=as.double(b[[5]]),as.double(con[[5]]),as.integer(cat[[5]]),sigma25=as.double(sigma2[[5]]),
    as.integer(c(ty[6],p[6],C[6])),alpha6=as.double(a[[6]]),beta6=as.double(b[[6]]),as.double(con[[6]]),as.integer(cat[[6]]),sigma26=as.double(sigma2[[6]]),
    as.integer(c(betaBurnin,betaDraw)),as.double(beta0),as.double(invSigma0),as.double(invSigmaBeta0),as.double(c(invga0,invgb0,betaVarScale)),as.double(pg),
    sumGamma1=as.integer(gamma[[1]]),sumABeta1=as.integer(acceptBeta[[1]]),sumAGa1=as.integer(acceptGamma[[1]]),
    sumGamma2=as.integer(gamma[[2]]),sumABeta2=as.integer(acceptBeta[[2]]),sumAGa2=as.integer(acceptGamma[[2]]),
    sumGamma3=as.integer(gamma[[3]]),sumABeta3=as.integer(acceptBeta[[3]]),sumAGa3=as.integer(acceptGamma[[3]]),
    sumGamma4=as.integer(gamma[[4]]),sumABeta4=as.integer(acceptBeta[[4]]),sumAGa4=as.integer(acceptGamma[[4]]),
    sumGamma5=as.integer(gamma[[5]]),sumABeta5=as.integer(acceptBeta[[5]]),sumAGa5=as.integer(acceptGamma[[5]]),
    sumGamma6=as.integer(gamma[[6]]),sumABeta6=as.integer(acceptBeta[[6]]),sumAGa6=as.integer(acceptGamma[[6]]),PACKAGE="iClusterPlus")

  Alpha = as.list(rep(NA,ndt))
  Beta =  as.list(rep(NA,ndt))
  Ratio =  as.list(rep(NA,ndt))
  acsGamma =  as.list(rep(NA,ndt))
  acsBeta =  as.list(rep(NA,ndt))
  Sigma2 = as.list(rep(NA,ndt))

  if(ndt >= 1){
    Alpha[[1]] = res$alpha1
    Beta[[1]] = matrix(res$beta1,ncol=K)
    Sigma2[[1]] = res$sigma21
    Ratio[[1]] = res$sumGamma1/betaDraw
    acsGamma[[1]] = res$sumAGa1/betaDraw
    acsBeta[[1]] = res$sumABeta1/betaDraw
    if(ty[1] == 1){
      acsBeta[[1]] = rep(1,p[1]) #normal data, acceptance ratio is always 1
    }
  }

  if(ndt >= 2){
    Alpha[[2]] = res$alpha2
    Beta[[2]] = matrix(res$beta2,ncol=K)
    Sigma2[[2]] = res$sigma22
    Ratio[[2]] = res$sumGamma2/betaDraw
    acsGamma[[2]] = res$sumAGa2/betaDraw
    acsBeta[[2]] = res$sumABeta2/betaDraw
    if(ty[2] == 1){
      acsBeta[[2]] = rep(1,p[2]) #normal data, acceptance ratio is always 1
    }  
  }

  if(ndt >= 3){
    Alpha[[3]] = res$alpha3
    Beta[[3]] = matrix(res$beta3,ncol=K)
    Sigma2[[3]] = res$sigma23
    Ratio[[3]] = res$sumGamma3/betaDraw
    acsGamma[[3]] = res$sumAGa3/betaDraw
    acsBeta[[3]] = res$sumABeta3/betaDraw
    if(ty[3] == 1){
      acsBeta[[3]] = rep(1,p[3]) #normal data, acceptance ratio is always 1
    }  
  }

  if(ndt >= 4){
    Alpha[[4]] = res$alpha4
    Beta[[4]] = matrix(res$beta4,ncol=K)
    Sigma2[[4]] = res$sigma24
    Ratio[[4]] = res$sumGamma4/betaDraw
    acsGamma[[4]] = res$sumAGa4/betaDraw
    acsBeta[[4]] = res$sumABeta4/betaDraw
    if(ty[4] == 1){
      acsBeta[[4]] = rep(1,p[4]) #normal data, acceptance ratio is always 1
    }      
  }

  if(ndt >= 5){
    Alpha[[5]] = res$alpha5
    Beta[[5]] = matrix(res$beta5,ncol=K)
    Sigma2[[5]] = res$sigma25
    Ratio[[5]] = res$sumGamma5/betaDraw
    acsGamma[[5]] = res$sumAGa5/betaDraw
    acsBeta[[5]] = res$sumABeta5/betaDraw
    if(ty[5] == 1){
      acsBeta[[5]] = rep(1,p[5]) #normal data, acceptance ratio is always 1
    }      
  }

  if(ndt >= 6){
    Alpha[[6]] = res$alpha6
    Beta[[6]] = matrix(res$beta6,ncol=K)
    Sigma2[[6]] = res$sigma26
    Ratio[[6]] = res$sumGamma6/betaDraw
    acsGamma[[6]] = res$sumAGa6/betaDraw
    acsBeta[[6]] = res$sumABeta6/betaDraw
    if(ty[6] == 1){
      acsBeta[[6]] = rep(1,p[6]) #normal data, acceptance ratio is always 1
    }      
  }
  
  list(meanZ=matrix(res$sumMeanZ,nrow=n,ncol=K),lastZ=matrix(res$lastZ,nrow=n,ncol=K),acsZ=res$acceptZ/zDraw/betaDraw,
       Alpha=Alpha, Beta=Beta, Sigma2=Sigma2,Ratio=Ratio,acsGamma=acsGamma,acsBeta=acsBeta)
}


iClusterBayes <- function(dt1,dt2=NULL,dt3=NULL,dt4=NULL,dt5=NULL,dt6=NULL,type = c("gaussian","binomial","poisson"),K=2,
                          n.burnin=1000,n.draw=1200,prior.gamma=rep(0.1,6),sdev=0.5,beta.var.scale=1,thin=1,pp.cutoff=0.5){

  dttype = c("gaussian","binomial","poisson")
  if(missing(dt1) | is.null(dt1)){
    stop("Error: dt1 is missing!\n")
  }

  if(!all(type %in% dttype)){
      cat("Error: ",type[!all(type %in% dttype)],"\n")
      stop("Allowed data types are gaussian, binomial and poisson. \n")
  }

  isNULL = c(is.null(dt1),is.null(dt2),is.null(dt3),is.null(dt4),is.null(dt5),is.null(dt6))
  if(any(diff(isNULL) == -1)){
      stop("Error: dt1 to dt6 must be assigned in order.\n")
  }

 if(sum(!isNULL) > length(type)){
     stop("Error:  data type is missing for some data. \n")
 } 
         
  ### burnin and draw for latent variable Z ##
  zBurnin = 0
  zDraw = 1
  ### the above settings make the draws for Z are the same as the draws for beta
  ### if zDraw > 1, the outcome Z from mcmcMix will be the mean Z of zDraw
  
  ### n.burnin and n.draw control the overall draw for Z and parameter beta
  betaBurnin = n.burnin
  betaDraw = n.draw
  
  n = nrow(dt1)
  ndt = 1
  Data = list()
  
  Data[[1]] = dataType(dt1,type[1],K)
   
  Data[[2]] = NULL
  if (!is.null(dt2)) {
    if(n != nrow(dt2)){
      stop("Error: nrow(dt1) != nrow(dt2) \n")
    }    
    Data[[2]] = dataType(dt2, type[2], K)
    ndt = ndt + 1
  }
  
  Data[[3]] = NULL
  if(!is.null(dt3)){
    if(n != nrow(dt3)){stop("Error: nrow(dt1) != nrow(dt3) \n")} 
    Data[[3]] = dataType(dt3,type[3],K)
    ndt = ndt + 1
  }
  
  Data[[4]] = NULL
  if(!is.null(dt4)){
    if(n != nrow(dt4)){stop("Error: nrow(dt1) != nrow(dt4) \n")} 
    Data[[4]] = dataType(dt4,type[4],K)
    ndt = ndt + 1
  }

  Data[[5]] = NULL
  if(!is.null(dt5)){
    if(n != nrow(dt5)){stop("Error: nrow(dt1) != nrow(dt5) \n")} 
    Data[[5]] = dataType(dt5,type[5],K)
    ndt = ndt + 1
  }

  Data[[6]] = NULL
  if(!is.null(dt6)){
      if(n != nrow(dt6)){stop("Error: nrow(dt1) != nrow(dt6) \n")} 
      Data[[6]] = dataType(dt6,type[6],K)
      ndt = ndt + 1
  }
  ###### priors for Bayesian variable selection #######
  invSigma0 = diag(rep(1,K+1))
  beta0 = rep(0,K+1)
  invSigmaBeta0 = invSigma0 %*% beta0
  invga0=1
  invgb0=1

  res = mcmcBayes(dt1=Data[[1]],dt2=Data[[2]],dt3=Data[[3]],dt4=Data[[4]],dt5=Data[[5]],dt6=Data[[6]],ndt,sdev,n,K,zBurnin,zDraw,
    betaBurnin, betaDraw,thin,prior.gamma,beta0,invSigma0,invSigmaBeta0,invga0,invgb0,beta.var.scale)

  for(i in 1:ndt){
    Data[[i]]$Alpha = res$Alpha[[i]]
    Data[[i]]$Beta = res$Beta[[i]]
    Data[[i]]$sigma2 = res$Sigma2[[i]]
    Data[[i]]$Ratio = res$Ratio[[i]]
    Data[[i]]$acsGamma = res$acsGamma
    Data[[i]]$acsBeta = res$acsBeta
  }
  BIC = totalBICbayes(Data, res$meanZ, ndt, K, pp.cutoff)
  devRatio = dev.ratio.bayes(Data, res$meanZ,ndt, K, pp.cutoff)
  kmeans.fit = kmeans(res$meanZ, K + 1, nstart=100)
  clusters = kmeans.fit$cluster
  centers = kmeans.fit$centers
  
  list(alpha=res$Alpha, beta=res$Beta, beta.pp=res$Ratio,gamma.ar=res$acsGamma,beta.ar=res$acsBeta,Z.ar=res$acsZ,
       clusters = clusters, centers = centers, meanZ = res$meanZ, BIC = BIC, dev.ratio = devRatio)
}


tune.iClusterBayes = function(cpus=6,dt1,dt2=NULL,dt3=NULL,dt4=NULL,dt5=NULL,dt6=NULL,type=c("gaussian","binomial","poisson"),
    K=1:6,n.burnin=1000,n.draw=1200,prior.gamma=rep(0.1,6),sdev=0.5,beta.var.scale=1,thin=1,pp.cutoff=0.5){

  #require(parallel) 
  dttype = c("gaussian","binomial","poisson")
  if(missing(dt1) | is.null(dt1)){
    stop("Error: dt1 is missing!\n")
  }

  if(!all(type %in% dttype)){
      cat("Error: ",type[!all(type %in% dttype)],"\n")
      stop("Allowed data types are gaussian, binomial and poisson. \n")
  }

  isNULL = c(is.null(dt1),is.null(dt2),is.null(dt3),is.null(dt4),is.null(dt5),is.null(dt6))
  if(any(diff(isNULL) == -1)){
      stop("Error: dt1 to dt6 must be assigned in order.\n")
  }

 if(sum(!isNULL) > length(type)){
     stop("Error:  data type is missing for some data. \n")
 } 
          
  cat("Begin parallel computation\n")
  RNGkind("L'Ecuyer-CMRG")
  fit = mclapply(K,FUN=function(x)iClusterBayes(dt1,dt2,dt3,dt4,dt5,dt6,type,x,n.burnin,n.draw,prior.gamma,sdev,beta.var.scale,thin,pp.cutoff), 
       mc.silent=TRUE, mc.cores=cpus, mc.preschedule=FALSE)
  cat("End parallel computation\n")	

  list(fit=fit)
}


plotHMBayes = function(fit, datasets, type = c("gaussian", "binomial", "poisson"),
    sample.order = NULL, row.order = NULL, sparse = NULL, 
    threshold = rep(0.5,length(datasets)), width = 5, scale = rep("none",length(datasets)), 
    col.scheme = rep(list(bluered(256)),length(datasets)),chr=NULL, plot.chr=NULL, cap=NULL)
{
    m = length(datasets)
    if(m > length(type)){
        stop("Error:  data type is missing for some data. \n")        
    }
    
    dttype = c("gaussian","binomial","poisson")
    if(!all(type %in% dttype)){
        cat("Error: ",type[!all(type %in% dttype)],"\n")
        stop("Allowed data types are gaussian, binomial and poisson. \n")
    }

    if (is.null(row.order)) {
        row.order = rep(T, m)
    }
    if (is.null(scale)) {
        scale = rep("none", m)
    }
    if (is.null(sparse)) {
        sparse = rep(F, m)
    }
    if (is.null(cap)) {
        cap = rep(F, m)
    }
    if (is.null(plot.chr)) {
        plot.chr = rep(F, m)
    }
    clusters = fit$clusters
    k = length(unique(clusters))
    if (is.null(sample.order)) {
        sorder = order(clusters)
    }
    else {
        sorder = sample.order
    }
    m = length(datasets)
    pp = unlist(lapply(1:m, function(l) {
        dim(datasets[[l]])[2]
    }))
    n = dim(datasets[[1]])[1]
    a = clusters[sorder]
    l = length(a)
    brkpoints = which(a[2:l] != a[1:(l - 1)])
    cluster.start = c(1, brkpoints + 1)
    my.panel.levelplot <- function(...) {
        panel.levelplot(...)
        panel.abline(v = (cluster.start[-1] - 0.5), col = "black", 
            lwd = 1, lty = 1)
        panel.scales = list(draw = FALSE)
    }
    for (i in 1:m) {
        pp = fit$beta.pp[[i]]
        upper = threshold[i]
        cat(i," ", sum(pp > upper),"\n")
        if (sparse[i] == T & sum(pp > upper) > 1) {
            image.data = datasets[[i]][sorder, which(pp > upper)]
        }else{
             warning("No variable selected!")
              image.data = datasets[[i]][sorder, ]
        }
        if (row.order[i] == T) {
            diss = 1 - cor(image.data, use = "na.or.complete")
            hclust.fit = hclust(as.dist(diss))
            gorder = hclust.fit$order
            image.data = image.data[, gorder]
        }
        if (plot.chr[i] == T) {
            if (sparse[i]) {
                chr = chr[which(pp > upper)]
            }
            len = length(chr)
            chrom.ends <- rep(NA, length(table(chr)))
            d = 1
            for (r in unique(chr)) {
                chrom.ends[d] <- max(which(chr == r))
                d = d + 1
            }
            chrom.starts <- c(1, chrom.ends[-length(table(chr))] + 
                1)
            chrom.mids <- (chrom.starts + chrom.ends)/2
            my.panel.levelplot.2 <- function(...) {
                panel.levelplot(...)
                panel.abline(v = (cluster.start[-1] - 0.5), col = "black", 
                  lwd = 1, lty = 1)
                panel.abline(h = len - chrom.starts[-1], col = "gray", 
                  lwd = 1)
                panel.scales = list(x = list(), y = list(at = len - 
                  chrom.mids), z = list())
            }
            my.panel = my.panel.levelplot.2
            scales = list(x = list(draw = F), y = list(at = len - 
                chrom.mids, labels = names(table(chr))), z = list(draw = F))
        }
        else {
            my.panel = my.panel.levelplot
            scales = list(draw = F)
        }
        scale.fn = function(x) {
            x <- sweep(x, 1L, rowMeans(x, na.rm = T), check.margin = T)
            sx <- apply(x, 1L, sd, na.rm = T)
            x <- sweep(x, 1L, sx, "/", check.margin = T)
            return(x)
        }
        if (scale[i] == "row") {
            image.data = scale.fn(image.data)
        }
        if (scale[i] == "col") {
            image.data = scale.fn(t(image.data))
            image.data = t(image.data)
        }
        image.data = as.matrix(rev(as.data.frame(image.data)))
        if (type[i] == "binomial") {
            colorkey = list(space = "right", height = 0.3, at = c(0, 
                0.5, 1), tick.number = 1)
        }
        else {
            colorkey = list(space = "right", height = 0.3, tick.number = 5)
        }
        if (cap[i] == T) {
            cut = quantile(datasets[[i]], prob = 0.9995, na.rm = T)
            p = levelplot(image.data, panel = my.panel, scales = scales, 
                col.regions = col.scheme[[i]], at = c(-Inf, seq(-cut, 
                  cut, length = 256), Inf), xlab = "", ylab = "", 
                colorkey = colorkey)
        }
        else {
            p = levelplot(image.data, panel = my.panel, scales = scales, 
                col.regions = col.scheme[[i]], xlab = "", ylab = "", 
                colorkey = colorkey)
        }
        if (i == m) {
            print(p, split = c(1, i, 1, m), more = F, panel.width = list(width, 
                "inches"))
        }
        else {
            print(p, split = c(1, i, 1, m), more = T, panel.width = list(width, 
                "inches"))
        }
    }
}

########### deviance ratio function using glmnet functions ############
dev.ratio.bayes0 = function(Data,meanZ,ndt,K){
  cat("- deviance ratio -\n")
  alpha = rep(1,ndt)
  lambda = rep(0,ndt) #set all lambda to zero because no-significant variables will be removed
  sumdev0 = 0
  sumdev = 0
  for(i in 1:ndt){
    #sigID = which(Data[[i]]$Ratio >= quantile(Data[[i]]$Ratio,prob=0.9))
    sigID = which(Data[[i]]$Ratio > 0.5)
    lenp = length(sigID)
    lenp0 = length(Data[[i]]$Ratio) - lenp
    if(Data[[i]]$type == 1){  # normal #
        if(lenp > 0){
            fit1 = .C("elnetBatchDev",a0=double(lenp),beta=double(lenp*K),sigma2 = double(lenp),
                as.double(meanZ),as.double(Data[[i]]$con[,sigID]),as.integer(Data[[i]]$n),as.integer(K),as.integer(lenp),
                as.double(alpha[i]),as.double(lambda[i]),sumdev0=double(1),sumdev=double(1))# ="iClusterPlus")
            cat(fit1$sumdev0, fit1$sumdev,"\n")
            sumdev0 = sumdev0 + fit1$sumdev0
            sumdev = sumdev + fit1$sumdev
        }
        if(lenp0 > 0){    
            fit1 = .C("elnetBatchDev",a0=double(lenp0),beta=double(lenp0*K),sigma2 = double(lenp0),
                as.double(meanZ),as.double(Data[[i]]$con[,-sigID]),as.integer(Data[[i]]$n),as.integer(K),as.integer(lenp0),
                as.double(alpha[i]),as.double(1),sumdev0=double(1),sumdev=double(1))# ="iClusterPlus")
            cat(fit1$sumdev0, fit1$sumdev,"\n")
            sumdev0 = sumdev0 + fit1$sumdev0
            sumdev = sumdev + fit1$sumdev
        }
    }else if(Data[[i]]$type == 2){ # binomial #
        if(lenp > 0){   
            fit2 = .C("lognetBatchDev",a0=double(lenp),beta=double(lenp*K),as.double(meanZ),
                as.integer(Data[[i]]$cat[,sigID]),as.integer(Data[[i]]$n),as.integer(K),as.integer(lenp),
                as.double(alpha[i]),as.double(lambda[i]),as.integer(Data[[i]]$nclass),as.integer(2),
                as.integer(0),sumdev0=double(1),sumdev=double(1))# ="iClusterPlus") #family=0 is binomial
            cat(fit2$sumdev0, fit2$sumdev,"\n")
            sumdev0 = sumdev0 + fit2$sumdev0
            sumdev = sumdev + fit2$sumdev
        }
        if(lenp0 > 0){
            fit2 = .C("lognetBatchDev",a0=double(lenp0),beta=double(lenp0*K),as.double(meanZ),
                as.integer(Data[[i]]$cat[,-sigID]),as.integer(Data[[i]]$n),as.integer(K),as.integer(lenp0),
                as.double(alpha[i]),as.double(1),as.integer(Data[[i]]$nclass),as.integer(2),
                as.integer(0),sumdev0=double(1),sumdev=double(1))# ="iClusterPlus") #family=0 is binomial
            cat(fit2$sumdev0, fit2$sumdev,"\n")
            sumdev0 = sumdev0 + fit2$sumdev0
            sumdev = sumdev + fit2$sumdev
        }
    }else if(Data[[i]]$type == 3){ # Poisson #
        if(lenp > 0){        
            fit3 = .C("fishnetBatchDev",a0=double(lenp),beta=double(lenp*K),as.double(meanZ),
                as.double(Data[[i]]$cat[,sigID]),as.integer(Data[[i]]$n),as.integer(K),as.integer(lenp),
                as.double(alpha[i]),as.double(lambda[i]),sumdev0=double(1),sumdev=double(1))# ="iClusterPlus")
            cat(fit3$sumdev0, fit3$sumdev,"\n")
            sumdev0 = sumdev0 + fit3$sumdev0
            sumdev = sumdev + fit3$sumdev
        }
        if(lenp0 > 0){
            fit3 = .C("fishnetBatchDev",a0=double(lenp0),beta=double(lenp0*K),as.double(meanZ),
                as.double(Data[[i]]$cat[,-sigID]),as.integer(Data[[i]]$n),as.integer(K),as.integer(lenp0),
                as.double(alpha[i]),as.double(1),sumdev0=double(1),sumdev=double(1))# ="iClusterPlus")
            cat(fit3$sumdev0, fit3$sumdev,"\n")
            sumdev0 = sumdev0 + fit3$sumdev0
            sumdev = sumdev + fit3$sumdev
        }
      }else { # Multinomial #
       # fit4 = .C("lognetBatchDev",a0=double(lenp * Data[[i]]$C),beta=double(lenp*K*Data[[i]]$C),
       #   as.double(meanZ),as.integer(Data[[i]]$cat[,sigID]),as.integer(Data[[i]]$n),as.integer(K),
       #   as.integer(lenp),as.double(alpha[i]),as.double(lambda[i]),as.integer(Data[[i]]$nclass),
       #   as.integer(Data[[i]]$C),as.integer(1),sumdev0=double(1),sumdev=double(1))# ="iClusterPlus") #family=1 is multinomial
       # sumdev0 = sumdev0 + fit4$sumdev0
       # sumdev = sumdev + fit4$sumdev
      }
  }
  # deviance ratio, for linear reg, it is R-square; the bigger, the better
  1-sumdev/sumdev0
}
