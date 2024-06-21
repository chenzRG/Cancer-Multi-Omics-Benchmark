###proportion of deviation from perfect block diagonal matrix###
compute.pod=function(fit){
  #deviation from block structure
  expZ=fit$meanZ
  cl=fit$cluster
  sorted=sort(table(cl))
  o.stol=as.numeric(names(sorted))
  o=NULL
  for(i in o.stol){
    o=c(o,which(cl==i))
  }
  a.matrix=t(expZ)%*%expZ
  diag.elements=diag(a.matrix)
  n=length(diag.elements)
  denom=matrix(rep(diag.elements,n),nrow=n, byrow=T)
  rr.matrix=a.matrix/sqrt(denom)/sqrt(t(denom))

  rr.matrix=replace(rr.matrix,rr.matrix<0,0)
  rr.matrix=rr.matrix[o,o]

  Z=model.matrix(~0+as.factor(cl))
  B=Z%*%t(Z)
  B=B[o,o]

  dev=sum(abs(B-rr.matrix))
  n=dim(expZ)[2]
  pod=dev/n/n
  return(pod)
}
