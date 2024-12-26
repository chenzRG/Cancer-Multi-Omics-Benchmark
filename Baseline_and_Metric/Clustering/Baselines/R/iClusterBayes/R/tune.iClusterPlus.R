#valid npoint value for the uniform design
#valid npoint for 2 lambdas are 8, 13, 21, 34, 55, 89, 144, 233, 377, 610
#valid npoint for 3 lambdas are 35, 101, 135, 185, 266, 418, 579, 828, 1010
#valid npoint for 4 lambdas are 307, 526, 701, 1019, 2129, 3001, 4001, 5003, 6007
#valid npoint for 5 lambdas are 1069, 1543, 2129, 3001, 4001, 5003, 6007, 8191
#no restriction for 1 lambda

#alpha and lambda are penalty parameter. alpha=1 correspons to lasso 0>alpha<1 elasticnet

tune.iClusterPlus = function(cpus=8,dt1,dt2=NULL,dt3=NULL,dt4=NULL,type=c("gaussian","binomial","poisson","multinomial"),
  K=2,alpha=c(1,1,1,1),n.lambda=NULL,scale.lambda=c(1,1,1,1),n.burnin=200,n.draw=200,maxiter=20,sdev=0.05,eps=1.0e-4){

#  require(parallel)
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
  
  lbd = list()
  lbd[[1]] = 1:1000
  lbd[[2]] = c(8, 13, 21, 34, 55, 89, 144, 233, 377, 610)
  lbd[[3]] = c(35, 101, 135, 185, 266, 418, 579, 828, 1010)
  lbd[[4]] = c(307, 526, 701, 1019, 2129, 3001, 4001, 5003, 6007)
  lbd[[5]] = c(1069,1543, 2129, 3001, 4001, 5003, 6007, 8191)

  ########### prepare sample points for lambda #################
  #### load("/home/nfs/moq/iClusterPlus/data/glp.rda") ###

  data(glp)
  m = 1
  m = m + as.numeric(!is.null(dt2)) + as.numeric(!is.null(dt3)) + as.numeric(!is.null(dt4))
  n.glp=c(50,144,185,307,1069)

  if(is.null(n.lambda)){
    n.lambda=n.glp[m]
    cat(n.glp[m]," points of lambdas are used to tune parameters.\n")
  }else{
    if(any(n.lambda==lbd[[m]])){
      cat(n.lambda," points of lambdas are used to tune parameters.\n")
    }else if(m>1){
      cat("Error: n.lambda is not a valid value\n")
      cat("Valid value for n.lambda are ",lbd[[m]],"\n")
      stop()
    }else{
      warning()
      cat(n.glp[m]," points of lambdas are used to tune parameters; may be too many.\n")
    }
  }
    
  which.glp=paste("s",m, ".txt", sep='')
  h=glp[[which.glp]]
  h=h[-1,which(h[1,]==n.lambda)]
  h=c(1,h)
  if(m==1){ud=as.matrix(seq(1,2*n.lambda-1, by=2)/2/n.lambda,nrow=n.lambda)}
  if(m>1){
    ud=matrix(NA,nrow=n.lambda, ncol=m)
    for(s in 1:m){
      ud[,s]=((1:n.lambda)*h[s]-0.5)/n.lambda	
    }
  }
  ud=ud-floor(ud)
  
  i=1
  while((i<=m) && (i<=length(scale.lambda))){
    ud[,i] = ud[,i]*scale.lambda[i]
    i = i+1
  }
     
   cat("Begin parallel computation\n")
   RNGkind("L'Ecuyer-CMRG")
   fit = mclapply(1:nrow(ud), 
       FUN=function(x)iClusterPlus(dt1, dt2,dt3,dt4,lambda=ud[x,],type, K, alpha,n.burnin,n.draw,maxiter, sdev, eps), 
       mc.silent=TRUE, mc.cores=cpus, mc.preschedule=FALSE)
  cat("End parallel computation\n")	

  list(fit=fit,lambda=ud)
}

