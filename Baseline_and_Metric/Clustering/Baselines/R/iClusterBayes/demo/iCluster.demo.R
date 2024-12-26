#setwd("/home/nfs/shenr/iClusterPLoSONEsoftware/")
library(iCluster)

data(breast.chr17)
fit=iCluster(datasets=breast.chr17, k=4, lambda=c(0.2,0.2))

plotiCluster(fit)

#a simple data transformation for nice image contrast
cn=breast.chr17[[2]]
cn[cn< -1.5]= -1.5;cn[cn>1.5]=1.5
exp=breast.chr17[[1]]
exp[exp< -2.5]= -2.5; exp[exp>2.5]=2.5

library(gplots)
library(lattice)
plotHeatmap(datasets=list(cn,exp),fit=fit)

#GBM data example including three data matrices: copy number, methylation, and expression.
data(gbm) 
data(coord)
chr=coord[,1]

#A variant of iCluster method with variance weighted shrinkage
fit=iCluster2(datasets=gbm, k=3, lambda=list(0.44,0.33,0.28))
             
plotiCluster(fit=fit)

compute.pod(fit)

plotHeatmap(fit=fit, datasets=gbm, feature.order=c(FALSE,TRUE,TRUE), sparse=c(FALSE,TRUE,TRUE),plot.chr=c(TRUE,FALSE,FALSE), chr=chr)

#Model tuning (using a simulated data example)
data(simu.datasets)
cv.fit=alist()
for(k in 2:5){
  cat(paste("K=",k,sep=""),'\n')
  cv.fit[[k]]=tune.iCluster2(datasets=simu.datasets, k,nrep=2, n.lambda=8)
}

#Reproducibility index (RI) plot
plotRI(cv.fit)


#Based on the RI plot, k=3 is the best solution
best.fit=cv.fit[[3]]$best.fit
#Try different color schemes
plotHeatmap(fit=best.fit,datasets=simu.datasets,sparse=c(TRUE,TRUE),col.scheme=list(bluered(256), greenred(256)))

