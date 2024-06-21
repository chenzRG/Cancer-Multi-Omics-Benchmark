
#library(lattice) 
#library(gplots)

#################################################
plotHeatmap=function(fit, datasets, type=c("gaussian","binomial","poisson","multinomial"), sample.order=NULL, row.order=NULL, sparse=NULL, threshold=rep(0.25,length(datasets)), width=5, scale=rep("none",length(datasets)), col.scheme=rep(list(bluered(256)),length(datasets)), chr=NULL, plot.chr=NULL, cap=NULL){

m=length(datasets)  
if(is.null(row.order)){row.order=rep(T,m)}  
if(is.null(scale)){scale=rep("none",m)}
if(is.null(sparse)){sparse=rep(F,m)}
if(is.null(cap)){cap=rep(F,m)}
if(is.null(plot.chr)){plot.chr=rep(F,m)}
  
#get clusters 
clusters=fit$clusters

k=length(unique(clusters))
if(is.null(sample.order)){sorder=order(clusters)}else{sorder=sample.order}
m=length(datasets)
pp=unlist(lapply(1:m,function(l){dim(datasets[[l]])[2]}))
n=dim(datasets[[1]])[1]


#cluster divider
a=clusters[sorder]
l=length(a)
brkpoints=which(a[2:l]!=a[1:(l-1)])
cluster.start=c(1,brkpoints+1)

my.panel.levelplot <- function(...) {
   panel.levelplot(...)
   panel.abline(v=(cluster.start[-1]-0.5), col="black", lwd=1,lty=1)
   panel.scales=list(draw=FALSE)
}



for(i in 1:m){
  
rowsum=apply(abs(fit$beta[[i]]),1, sum)
if(sum(rowsum)==0)warning(paste("All Lasso coefficients are zero for data type",i), call. = FALSE)

if(mean(rowsum>0)>threshold[i]){upper=quantile(rowsum,prob=(1-threshold[i]))}else{upper=0}

if(sparse[i]==T&sum(rowsum> upper)>1){
	image.data=datasets[[i]][sorder,which(rowsum> upper)]
}else{image.data=datasets[[i]][sorder,]}  
  

if(row.order[i]==T){
  diss=1-cor(image.data,use="na.or.complete");
  hclust.fit=hclust(as.dist(diss));
  #hclust.fit=hclust(dist(t(image.data)))
  gorder=hclust.fit$order
  image.data=image.data[,gorder]
}

if(plot.chr[i]==T){
if(sparse[i]){chr=chr[which(rowsum > upper)]}
len=length(chr)
chrom.ends <- rep(NA,length(table(chr)));   
d=1
for(r in unique(chr))
  {
    chrom.ends[d] <-max(which(chr==r))
    d=d+1
  }
chrom.starts <-c(1,chrom.ends[-length(table(chr))]+1)
chrom.mids <- (chrom.starts+chrom.ends)/2


my.panel.levelplot.2 <- function(...) {
   panel.levelplot(...)
   panel.abline(v=(cluster.start[-1]-0.5), col="black", lwd=1,lty=1)
   panel.abline(h=len-chrom.starts[-1],col="gray", lwd=1)
   #labels=list(at=len-chrom.mids,cex=0.6,labels=names(table(chr)))
   panel.scales=list(x=list(), y=list(at=len-chrom.mids), z=list())
   
}

my.panel=my.panel.levelplot.2 

scales=list(x=list(draw=F), y=list(at=len-chrom.mids,labels=names(table(chr))), z=list(draw=F))}else{
my.panel=my.panel.levelplot; 
scales=list(draw=F)}


scale.fn=function(x){
x <- sweep(x, 1L, rowMeans(x, na.rm = T), check.margin = T)
sx <- apply(x, 1L, sd, na.rm = T)
x<- sweep(x, 1L, sx, "/", check.margin = T)
return(x)
}

if(scale[i]=='row'){image.data=scale.fn(image.data)}
if(scale[i]=='col'){image.data=scale.fn(t(image.data)); image.data=t(image.data)}

image.data=as.matrix(rev(as.data.frame(image.data)))   #reverse the rows so the image will be top-down ordering

if(type[i]=="binomial"){colorkey=list(space="right", height=0.3, at=c(0,0.5,1), tick.number=1)}else{colorkey=list(space="right", height=0.3, tick.number=5)}

if(cap[i]==T){cut=quantile(datasets[[i]],prob=0.9995,na.rm=T);
p=levelplot(image.data, panel = my.panel, scales=scales, col.regions = col.scheme[[i]],at = c(-Inf, seq(-cut, cut, length=256), Inf), xlab="",ylab="", colorkey=colorkey )}else
{p=levelplot(image.data, panel = my.panel, scales=scales, col.regions = col.scheme[[i]], xlab="",ylab="", colorkey=colorkey )}


#print(p,panel.width=list(4,"inches"))

if(i==m){
print(p,split=c(1,i,1,m),more=F,panel.width=list(width,"inches"))}else{
print(p,split=c(1,i,1,m),more=T, panel.width=list(width,"inches"))}

}

}
