###04/12: GenomicRanges findOverlaps slots changed from matchMatrix to queryHits and subjectHits  
#library(GenomicRanges)

########CNregions##############################################

#seg file format
# sample chromosome    start      end num.mark seg.mean
# 01         1         59839  9523923       1813       0.0111
# 02         1      10453531 10465787          6      -0.2691
# 03         1      12073884 12405989         77       1.7106
# 04         1      12412402 22280160       1467      -0.0188
# 05         1      22708838 23239214         43       0.2741
# 06         1      23254415 41017098       2454      -0.0165

#cnv file format
#    chr    start      end
#  chr1  1093942  1293942
#  chr1  1619653  1730263
#  chr1  8746272  8890240
#  chr1 13372572 13557162
#  chr1 17007592 17125658
#  chr1 37942158 38054381

#recommend epsilon=1/n

CNregions=function(seg, epsilon=0.005, adaptive=FALSE, rmCNV=FALSE, cnv=NULL, frac.overlap=0.5, rmSmallseg=TRUE, nProbes=15){
  
  colnames(seg)=c("sample","chromosome","start","end","num.mark","seg.mean")
  seg = subset(seg, chromosome<=22)
  seg = seg[order(seg[,1],seg[,2],seg[,3]), ]

  #this can be effective in removing cnvs that do not yet exist in the cnv database
  if(rmSmallseg){seg=subset(seg, num.mark>=nProbes)}

  ############
  #remove CNV#
  ############
  if(rmCNV){
      cat("Removing CNV...",'\n')
      if(is.null(cnv))stop("To remove CNV, please include the cnv file. Otherwise set rmCNV=F")
      gr.cnv = GRanges(seqnames=as.character(cnv[,1]),
             ranges=IRanges(start=as.numeric(cnv[,2]),
                            end=as.numeric(cnv[,3])))

      gr.seg = GRanges(seqnames=paste("chr",seg[,2],sep=""),
             ranges=IRanges(start=seg[,3],
                            end=seg[,4]))

      overlap=findOverlaps(gr.cnv, gr.seg)
      queryHits=queryHits(overlap)
      subjectHits=subjectHits(overlap)
      start=pmax(seg[subjectHits,3],cnv[queryHits,2])
      end=pmin(seg[subjectHits,4],cnv[queryHits,3])
      seg.length=seg[subjectHits,4]-seg[subjectHits,3]
      seg.length[seg.length==0]=1
      p.overlap=(end-start)/seg.length
      
      #seg: --sxxxxxxxxxxe--------
      #cnv: ---------+++++-------
      #overlap=5
      #p.overlap=5/10=0.5

      cnv.idx=subjectHits[which(p.overlap>= frac.overlap)]
      cnv.idx=unique(cnv.idx)

      seg=seg[-cnv.idx,]
  }

  ##########################
  #compile unique CNregions#
  ##########################
  if(dim(seg)[1]==0)stop("seg file is empty.")
  samples=unique(seg$sample)
  chr=start.maploc=end.maploc=nur=out=NULL
  for(i in unique(as.numeric(seg[,2]))){
  
      #cat('Processing chr',i,'\n')
  
      #get union of breakpoints from CBS seg file
      subdata=subset(seg, chromosome==i)
      #u.breakpts=c(sort(unique(subdata$start)),max(subdata$end))  
      u.breakpts=sort(unique(subdata[,3]))
      #sample1:s--------e-s--------e--s---e---s--e
      #sample2:s--e-s-----e-s--e--s-------e-s----e
      l=length(u.breakpts)
      outi=NULL
      for(j in 1:length(samples)){
          sj=subset(subdata,sample==samples[j]) 
          if(dim(sj)[1]==0)stop("Found sample(s) with no segments. Check parameter setting. Relax threshold values.")
          outj=rep(NA,l)   
          sj$start[1]=u.breakpts[1] #for samples starting position lags, force them to conform
          idx=c(match(sj$start,u.breakpts), l)
          outj=c(rep(sj[,6],times=diff(idx)), tail(sj[,6],1))  
          outi=cbind(outi, outj)
      }
    u.start=u.breakpts
    #s.end=sort(unique(subdata$end))
    u.se= sort(union(unique(subdata$end),unique(subdata$start))) 
    u.end=lapply(2:length(u.start),FUN=function(x)u.se[which(u.start[x]<=u.se)[1]])
    u.end=unlist(u.end)  
    #u.end=sort(unique(subdata$end))[unlist(u.end)]
    u.end=c(u.end, max(u.se))  
          
    #compute mean sqaured distance between adjacent rows
    adjrow.dist=apply(sqrt(diff(outi)^2),1,mean) 
    #add row 1-1 to recover the original length
    adjrow.dist=c(0, adjrow.dist) 
    if(adaptive){epsilon=quantile(adjrow.dist,prob=0.97)}
    seg.break=which(adjrow.dist>= epsilon)
    if(length(seg.break)==0)stop("Check parameter setting. Set lower epsilon value.")  
    start=seg.break
    if(start[1]>1){start=c(1,start)}
    end=seg.break-1
    if(end[1]==0){end=end[-1]}
    len=dim(outi)[1]
    if(tail(end,1)<len){end=c(end,len)}  
    nuri=end-start+1  #number of unique regions
    rowidx=rep(c(1:length(nuri)), nuri)
    nur=c(nur,nuri)
    
    #get the medoid signature of each region with RMSE< threshold 
    get.medoid=function(x){
      if(dim(x)[1]>1){pam(x,k=1)$medoid}else{x}
      }
    outi.medoid=by(outi,rowidx,get.medoid)
    outi.medoid=matrix(unlist(outi.medoid),nrow=length(outi.medoid),byrow=T)
          
    start.maploc=c(start.maploc,u.start[start])
    end.maploc=c(end.maploc, u.end[end])
    chr=c(chr,rep(i,dim(outi.medoid)[1]))
  	out=rbind(out,outi.medoid)
  
	}
	colnames(out)=samples
  rownames(out)=paste("chr", chr, ".",start.maploc,"-",end.maploc, sep="")
	reducedM=t(out)
  return(reducedM)
}

