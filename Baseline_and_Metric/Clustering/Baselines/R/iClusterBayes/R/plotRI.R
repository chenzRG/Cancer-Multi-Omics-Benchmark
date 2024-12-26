plotRI=function(cv.fit){  
  ri=unlist(lapply(2:5,FUN=function(x)max(cv.fit[[x]]$RI)))
  k=as.character(2:length(cv.fit))
  plot(k,ri, axes=F, ylab="Reproducibility index", xlab="K",type="b",pch=20)
  axis(side=1,at=k,labels=k)
  axis(side=2)
}