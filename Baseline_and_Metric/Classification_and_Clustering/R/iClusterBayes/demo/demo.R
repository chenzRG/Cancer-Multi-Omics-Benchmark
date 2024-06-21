#nohup R CMD BATCH --no-save --no-restore demo.R demo.Rout &

# code to generate simulated data sets
n1=20
n2=20
n3=20
n = n1+n2+n3
k=2
p=200
p1 = 10
p2 = 10
p3 = 10 

true.class=c(rep(1,n1),rep(2,n2),rep(3,n3))

set.seed(339)
################# normal distribution ####################
c1=matrix(rnorm(n1*p1,mean=3,sd=1), ncol=p1,nrow=n1)
c2=matrix(rnorm(n2*p2,mean=2, sd=1), ncol=p2,nrow=n2)
c3=matrix(rnorm(n3*p3,mean=1, sd=1), ncol=p3,nrow=n3)

normData = matrix(rnorm(n*p, mean=0, sd=1),nrow=n, ncol=p)
normData[1:n1,1:p1]= c1
normData[(n1+1):(n1+n2),(p1+1):(p1+p2)]= c2
normData[(n1+n2+1):n,(p1+p2+1):(p1+p2+p3)]= c3

######## binomial data #########################
c1=matrix(rbinom(n1*p1,size=1, prob=0.8), ncol=p1,nrow=n1)
c2=matrix(rbinom(n2*p2,size=1, prob=0.5), ncol=p2,nrow=n2)
c3=matrix(rbinom(n3*p3,size=1, prob=0.3), ncol=p3,nrow=n3)

binomData = matrix(rbinom(n*p, size=1, prob=0.05),nrow=n, ncol=p)
binomData[1:n1,1:p1]= c1
binomData[(n1+1):(n1+n2),(p1+1):(p1+p2)]= c2
binomData[(n1+n2+1):n,(p1+p2+1):(p1+p2+p3)]= c3

######### poisson data ############################
c1=matrix(rpois(n1*p1,lambda=5), ncol=p1,nrow=n1)
c2=matrix(rpois(n2*p2,lambda=3), ncol=p2,nrow=n2)
c3=matrix(rpois(n3*p3,lambda=1), ncol=p3,nrow=n3)

poisData = matrix(rpois(n*p, lambda=0.5),nrow=n, ncol=p)
poisData[1:n1,1:p1]= c1
poisData[(n1+1):(n1+n2),(p1+1):(p1+p2)]= c2
poisData[(n1+n2+1):n,(p1+p2+1):(p1+p2+p3)]= c3

################## multinomial distribution ##############
c1=matrix(sample(1:3,size=n1*p1, replace=TRUE, prob=c(0.8,0.15,0.05)), ncol=p1,nrow=n1)
c2=matrix(sample(1:3,size=n2*p2, replace=TRUE, prob=c(0.3,0.4,0.3)), ncol=p2,nrow=n2)
c3=matrix(sample(1:3,size=n3*p3, replace=TRUE, prob=c(0.05,0.15,0.8)), ncol=p3,nrow=n3)

multData = matrix(sample(1:3,size=n*p, replace=TRUE, prob=c(0.05,0.9,0.05)),nrow=n, ncol=p)
multData[1:n1,1:p1]= c1
multData[(n1+1):(n1+n2),(p1+1):(p1+p2)]= c2
multData[(n1+n2+1):n,(p1+p2+1):(p1+p2+p3)]= c3

#biID = seq(1,p,2)
#multData[,biID] = binomData[,biID]

### some variables only have one category; need to remove these variables
ncat = function(x){length(unique(x))}
len = apply(binomData,2,ncat)
unicat = which(len==1)
binomData = binomData[,-unicat]

simuResult = list()
for(i in 1:5){
  simuResult[[i]] = tune.iClusterPlus(cpus=8,clusterType="SOCK",dt1=normData,dt2=binomData,
              dt3=poisData,dt4=multData,type=c("gaussian","binomial","poisson","multinomial"),
              K=i,alpha=c(1,1,1,1),n.lambda=307,scale.lambda=c(0.5,0.5,0.5,0.5),
              n.burnin=200,n.draw=200,maxiter=20,sdev=0.05,eps=1.0e-4)
}

save(simuResult,file="simuResult.RData")

load("simuResult.RData")

nLambda = nrow(simuResult[[1]]$lambda)
nK = length(resultList)

resultList = simuResult

getBIC = function(resultList){
  n.lambda = nrow(resultList[[1]]$lambda)
  n.K = length(resultList)
  BIC = matrix(NA,nrow=n.lambda,ncol=n.K)
  for(i in 1:n.K){
    BIC[,i] = unlist(lapply(1:n.lambda, FUN=function(x){resultList[[i]]$fit[[x]]$BIC}))
  }
  colnames(BIC)=paste("K=",1:n.K,sep="")
  BIC
}

getDevR = function(resultList){
  n.lambda = nrow(resultList[[1]]$lambda)
  n.K = length(resultList)
  DEVR = matrix(NA,nrow=n.lambda,ncol=n.K)
  for(i in 1:n.K){
    DEVR[,i] = unlist(lapply(1:n.lambda, FUN=function(x){resultList[[i]]$fit[[x]]$dev.ratio}))
  }
  colnames(DEVR)=paste("K=",1:n.K,sep="")
  DEVR
}

BIC = getBIC(resultList)

devR = getDevR(resultList)

minBICid = apply(BIC,2,which.min)
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}

plot(1:nK,devRatMinBIC,type="b",xlab="K",ylab="Dev.ratio at minBIC")

getClusters = function(resultList){
  BIC = getBIC(resultList)
  minBICid = apply(BIC,2,which.min)
  clusters = matrix(NA,nrow=length(resultList[[1]]$fit[[1]]$clusters),ncol=ncol(BIC))
  for(i in 1:ncol(BIC)){
    clusters[,i] = resultList[[i]]$fit[[minBICid[i]]]$clusters    
  }
  clusters
}

getClusters(resultList)


par(mfrow=c(2,1))
plot(resultList[[2]]$fit[[minBICid[2]]]$beta[[1]][,1],xlab="p",ylab="value")
plot(resultList[[2]]$fit[[minBICid[2]]]$beta[[1]][,2],xlab="p",ylab="value")

         
pdf("devR.pdf",height=11.5,width=8.5)
plot(1:nK,devR[1,],ylim=c(min(devR,na.rm=TRUE),max(devR,na.rm=TRUE)),col=0,xlab="Lambda index",ylab="BIC",type="b")
for(i in 2:nrow(devR)){
  lines(1:nK,devR[i,],type="b",col=round(i/30,0))
}
dev.off()


plot(1:nK,BIC[1,],ylim=c(min(BIC,na.rm=TRUE),max(BIC,na.rm=TRUE)),col=0,xlab="Lambda index",ylab="BIC",type="b")
for(i in 2:nrow(BIC)){
  lines(1:nK,BIC[i,],type="b",col=round(i/30,0))
}


