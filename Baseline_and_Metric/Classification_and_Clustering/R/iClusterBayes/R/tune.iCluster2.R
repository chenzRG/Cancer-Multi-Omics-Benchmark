#validat n.lambda value for the uniform design
#valid n.lambda for 2 lambdas are 8, 13, 21, 34, 55, 89, 144, 233, 377, 610
#valid n.lambda for 3 lambdas are 35, 101, 135, 185, 266, 418, 579, 828, 1010
#valid n.lambda for 4 lambdas are 307, 526, 701, 1019, 2129, 3001, 4001, 5003, 6007
#valid n.lambda for 5 lambdas are 1069, 1543, 2129, 3001, 4001, 5003, 6007, 8191
#no restriction for 1 lambda

#revised 10/09/10 per sijian's suggestion on predict.kmeans
# last updated 10/09/2020
#library(mclust)

tune.iCluster2=function(x, K, method=c("lasso","enet","flasso","glasso","gflasso"), base=200, chr=NULL, true.class=NULL, lambda=NULL, n.lambda=NULL, save.nonsparse=F, nrep=10, eps=1e-4)
  {
    datasets=x
    n <-dim(datasets[[1]])[1]
    m=length(datasets) 
    p=unlist(lapply(1:m,function(l){dim(datasets[[l]])[2]}))
    sum.p=sum(p)

    #pre-checks
    if(is.list(datasets)==F)stop("datasets must be a list")
    if(sum(is.na(datasets))>0)stop("data cannot have NAs. please exclude or impute missing values.")
    if(missing(K))stop("must specify the number of clusters k")

    if(is.null(method)){method=rep("lasso",m)}
    nlam=ifelse(method%in%c("lasso","glasso"),1,ifelse(method%in%c("enet","flasso","gflasso"),2,NA))
    total.nlam=sum(nlam)
    if(is.na(total.nlam))stop("Method not correctly specified")

    if(is.null(lambda)){     
    lbd = list()
    lbd[[1]] = 1:1000
    lbd[[2]] = c(8, 13, 21, 34, 55, 89, 144, 233, 377, 610)
    lbd[[3]] = c(35, 101, 135, 185, 266, 418, 579, 828, 1010)
    lbd[[4]] = c(307, 526, 701, 1019, 2129, 3001, 4001, 5003, 
                 6007)
    lbd[[5]] = c(1069, 1543, 2129, 3001, 4001, 5003, 6007, 8191)
    data(glp)
    n.glp = c(50, 144, 185, 307, 1069)
    if (is.null(n.lambda)) {
      n.lambda = n.glp[m]
      cat(n.glp[m], " points of lambdas are used to tune parameters.\n")
    }else{
      if (any(n.lambda == lbd[[m]])) {
        cat(n.lambda, " points of lambdas are used to tune parameters.\n")
      }else if (m > 1) {
        cat("Error: n.lambda is not a valid value\n")
        cat("Valid value for n.lambda are ", lbd[[m]], "\n")
        stop()
      }else {
        warning()
        cat(n.glp[m], " points of lambdas are used to tune parameters; may be too many.\n")
      }
    }
    which.glp=paste("s",total.nlam, ".txt", sep='')
    h=glp[[which.glp]]
    h=h[-1,which(h[1,]==n.lambda)]
    h=c(1,h)
    if(total.nlam==1){ud=as.matrix(seq(1,2*n.lambda-1, by=2)/2/n.lambda,nrow=n.lambda)}
    if(total.nlam>1){
    	ud=matrix(NA,nrow=n.lambda, ncol=total.nlam)
        for(s in 1:total.nlam){
		ud[,s]=((1:n.lambda)*h[s]-0.5)/n.lambda
	}
    }
    ud=ud-floor(ud)
    ud=base^sqrt(ud)
    }else{
      ud=lambda
    }

   temp=0
   ps=ps.adjusted=pred.error=NULL
   for(c in 1:dim(ud)[1]){
       
    cat(paste("uniform sampling point",c),"\n")
    lambda.c=alist()
    a=cumsum(c(1,nlam))[-(m+1)]
    b=cumsum(nlam)
    for(l in 1:m){
	#if(method[l]%in%c("flasso","gflasso")){cc=ud[c,a[l]:b[l]]; lambda.c[[l]]=c(cc[1],log(cc[2],base)*10)}else #the smoothness lambda max set at 10
	{lambda.c[[l]]=ud[c,a[l]:b[l]]}	
    }

    quick.check=iCluster2(x=datasets, K=K, method=method, chr=chr, lambda=lambda.c, eps=eps)

    
    #if(sum(abs(apply(quick.check$Bmat,2,sum))<1e-8)>0){ps[c]=ps.adjusted[c]=pred.error[c]=NA}else{   #10/09/2020
    if(sum(abs(apply(do.call(rbind.data.frame, quick.check$beta),2,sum))<1e-8)>0){ps[c]=ps.adjusted[c]=pred.error[c]=NA}else{
       ps.c=ps.adjusted.c=pred.error.c=NULL    
    	for(nsplit in 1:nrep){
    	#cat(paste("data split", nsplit),"\n")
    	folds <- split(sample(seq(n)), rep(1:2, length=n))

      omit <- folds[[1]]
	    trdata=alist()
	    tsdata=alist()
	      for(j in 1:m){
        	trdata[[j]] <- datasets[[j]][-omit, ]
        	tsdata[[j]] <- datasets[[j]][omit, ]
        }
	    fit.trdata=iCluster2(x=trdata, K=K, method=method, chr=chr, lambda=lambda.c, eps=eps)
            #W=fit.trdata$Bmat                             #10/09/2020
            W = do.call(rbind.data.frame, fit.trdata$beta) #10/09/2020
            W=apply(W,2,as.numeric)                        #10/09/2020
	    PSI=fit.trdata$Phivec
	    sigma=W%*%t(W)+diag(PSI)
  	  stacked.tsdata=matrix(NA, nrow=length(omit), ncol=sum.p)
      	for (t in 1:m) {
     	     if(t==1){idx=1:p[t]}else{idx=(sum(p[1:(t-1)])+1):sum(p[1:t])}
     	     stacked.tsdata[,idx] <- tsdata[[t]]
      	 }
	    pred.expZ=t(W)%*%solve(sigma)%*%t(stacked.tsdata)


    	fit=iCluster2(x=tsdata,K=K,method=method, chr=chr, lambda=lambda.c, eps=eps)    
    	pred=predict.kmeans(fit.trdata,t(pred.expZ))
   	  ps.c[nsplit]=RandIndex(fit$cluster, pred)
   	  ps.adjusted.c[nsplit]=aRandIndex(fit$cluster, pred) 
    	if(!is.null(true.class)){pred.error.c[nsplit]=classError(true.class[omit],pred)$errorRate}
   	    
         }

  	  ps[c]=mean(ps.c)
  	  ps.adjusted[c]=mean(ps.adjusted.c)
          #pred.error[c]=mean(pred.error.c)
  	  if(!is.null(true.class)){pred.error[c]=mean(pred.error.c)}
	}

}

  if(sum(is.na(ps.adjusted))==length(ps.adjusted))stop("lambda out of range")
  best.ps.adjusted=max(ps.adjusted,na.rm=T)
  bl=ud[max(which(ps.adjusted>=(best.ps.adjusted-0.1))),]  #break ties by chossing the most sparse model
  best.lambda=alist()
    a=cumsum(c(1,nlam))[-(m+1)]
    b=cumsum(nlam)
    for(l in 1:m){
	{best.lambda[[l]]=bl[a[l]:b[l]]}	
    }


  best.fit=iCluster2(x=datasets,K=K, method=method, chr=chr, lambda=best.lambda,eps=eps)
  if(!is.null(true.class)){best.error=classError(true.class,best.fit$cluster)$errorRate}else{best.error=NULL}
  
  if(save.nonsparse){fit0=iCluster2(x=datasets,K=K, method=rep("lasso",m), lambda=rep(0,m),eps=eps)}else{fit0=NULL}
  out=list(best.fit=best.fit, best.lambda=best.lambda, best.error=best.error, pred.error=pred.error, nonsparse.fit=fit0, true.class=true.class, ud=ud, ps=ps, ps.adjusted=ps.adjusted)
  return(out)

}

aRandIndex=function (x, y) 
{
    x <- as.vector(x)
    y <- as.vector(y)
    xx <- outer(x, x, "==")
    yy <- outer(y, y, "==")
    upper <- row(xx) < col(xx)
    xx <- xx[upper]
    yy <- yy[upper]
    a <- sum(as.numeric(xx & yy))
    b <- sum(as.numeric(xx & !yy))
    c <- sum(as.numeric(!xx & yy))
    d <- sum(as.numeric(!xx & !yy))
    ni <- (b + a)
    nj <- (c + a)
    abcd <- a + b + c + d
    q <- (ni * nj)/abcd
    (a - q)/((ni + nj)/2 - q)
}

RandIndex=function (c1, c2)
{
    c1 <- as.vector(c1)
    c2 <- as.vector(c2)
    xx <- outer(c1, c1, "==")
    yy <- outer(c2, c2, "==")
    upper <- row(xx) < col(xx)
    xx <- xx[upper]
    yy <- yy[upper]
    a <- sum(as.numeric(xx & yy))
    d <- sum(as.numeric(!xx & !yy))
    (a+d)/choose(length(c2),2)
}

classError=function (classification, truth) 
{
  q <- function(map, len, x) {
    x <- as.character(x)
    map <- lapply(map, as.character)
    y <- sapply(map, function(x) x[1])
    best <- y != x
    if (all(len) == 1) 
      return(best)
    errmin <- sum(as.numeric(best))
    z <- sapply(map, function(x) x[length(x)])
    mask <- len != 1
    counter <- rep(0, length(len))
    k <- sum(as.numeric(mask))
    j <- 0
    while (y != z) {
      i <- k - j
      m <- mask[i]
      counter[m] <- (counter[m]%%len[m]) + 1
      y[x == names(map)[m]] <- map[[m]][counter[m]]
      temp <- y != x
      err <- sum(as.numeric(temp))
      if (err < errmin) {
        errmin <- err
        best <- temp
      }
      j <- (j + 1)%%k
    }
    best
  }
  if (any(isNA <- is.na(classification))) {
    classification <- as.character(classification)
    nachar <- paste(unique(classification[!isNA]), collapse = "")
    classification[isNA] <- nachar
  }
  MAP <- mapClass(classification, truth)
  len <- sapply(MAP[[1]], length)
  if (all(len) == 1) {
    CtoT <- unlist(MAP[[1]])
    I <- match(as.character(classification), names(CtoT), 
               nomatch = 0)
    one <- CtoT[I] != truth
  }
  else {
    one <- q(MAP[[1]], len, truth)
  }
  len <- sapply(MAP[[2]], length)
  if (all(len) == 1) {
    TtoC <- unlist(MAP[[2]])
    I <- match(as.character(truth), names(TtoC), nomatch = 0)
    two <- TtoC[I] != classification
  }
  else {
    two <- q(MAP[[2]], len, classification)
  }
  err <- if (sum(as.numeric(one)) > sum(as.numeric(two))) 
    as.vector(one)
  else as.vector(two)
  bad <- seq(along = classification)[err]
  list(misclassified = bad, errorRate = length(bad)/length(truth))
}

predict.kmeans=function(km, data)
{k <- nrow(km$centers)
n <- nrow(data)
d <- as.matrix(dist(rbind(km$centers, data)))[-(1:k),1:k]
out <- apply(d, 1, which.min)
return(out)}
