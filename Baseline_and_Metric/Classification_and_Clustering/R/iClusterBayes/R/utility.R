
# get the BIC for each lambda and K; rows correspond to the lambda (vector)
# and columns correspond to the K latent variables

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

# get the deviance ratio for each lambda and K; rows correspond to the lambda (vector)
# and columns correspond to the K latent variables
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

### get the cluster assignments for 
getClusters = function(resultList){
  BIC = getBIC(resultList)
  minBICid = apply(BIC,2,which.min)
  clusters = matrix(NA,nrow=length(resultList[[1]]$fit[[1]]$clusters),ncol=ncol(BIC))
  for(i in 1:ncol(BIC)){
    clusters[,i] = resultList[[i]]$fit[[minBICid[i]]]$clusters    
  }
  clusters
}

iManual = function (view = TRUE) 
{
    ipdf <- system.file("doc", "iManual.pdf", package = "iClusterPlus")
    if (view) {
        if (.Platform$OS.type == "windows") 
            shell.exec(ipdf)
        else system(paste(Sys.getenv("R_PDFVIEWER"),ipdf, "&"))
    }
    return(ipdf)
}
