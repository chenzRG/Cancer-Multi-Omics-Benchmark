plotiCluster = function(fit,label=NULL){
  
   cl=fit$clusters
   #order clusters from smallest cluster to largest
   sorted=sort(table(cl))
   o.stol=as.numeric(names(sorted))
   o=NULL
   for(i in o.stol){
   o=c(o,which(cl==i))
   }
   
   s.matrix=t(fit$meanZ)%*%fit$meanZ
   #standardize s.matrix such that the diagonal is 1
   diag.elements=diag(s.matrix)
   n=length(diag.elements)
   denom=matrix(rep(diag.elements,n),nrow=n, byrow=T)
   a=s.matrix/sqrt(denom)/sqrt(t(denom))
   a=replace(a,a<0,0) #assuming negative correlation is not meaningful

   a=a[o,o]
   n=dim(a)[1]
   
   #flip matrix for plot orientation
   f.a=t(as.matrix(rev(as.data.frame(t(a)))))

   image(1:(ncol(f.a)+1), 1:(nrow(f.a)+1),t(f.a), axes=FALSE, col=gray(25:0/25),ylab="",xlab="")
   if(!is.null(label)){
     axis(side=1, at=1.5:(n+0.5), label[o],las=2,cex.axis=0.5)
     axis(side=2, at=(n+0.5):1.5, label[o],las=2,cex.axis=0.5)
     }else{
     axis(side=1, at=1:n, 1:n)
     axis(side=2, at=1:n, 1:n)
     }
   mtext(side=3,line=1,paste("K=",max(cl),sep=""))
   box()

}
