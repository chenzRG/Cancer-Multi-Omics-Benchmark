setwd('D:/brca')
library('plyr')
library("ipred")
library("survival")
library("survivalROC")
library("glmnet")
kx=read.csv("brca_cox.csv",row.names= 1)
mysurv=Surv(kx$time,kx$status)


Unicox <- function(x){
  #cat("aa:",x,"\n")
  fml <- as.formula(paste0('mysurv~', x))
  gcox <- coxph(fml, kx)
  cox_sum <- summary(gcox)
  HR <- round(cox_sum$coefficients[,2],2)
  PValue <- round(cox_sum$coefficients[,5],4)
  CI <- paste0(round(cox_sum$conf.int[,3:4],2),collapse='-')
  Uni_cox <- data.frame('Characteristics' = x,
                        'Hazard Ratio' = HR,
                        'CI95' = CI,
                        'P value' = PValue)
  return(Uni_cox)
}
VarNames <-colnames(kx)[3:ncol(kx)]
Univar <- lapply(VarNames, Unicox)
Univar <- ldply(Univar, data.frame)

Univar[,5]=p.adjust(Univar$P.value, method ="fdr", n=dim(Univar)[1])#
colnames(Univar)<-c("Characteristics", "Hazard.Ratio", "CI95", "P.value", "adj.p")
dd=matrix(0,1,ncol(kx))
dd[1,3:ncol(kx)]=Univar$P.value
dd=data.frame(dd)
colnames(dd)=colnames(kx)
ee=rbind(dd,kx)
ff=ee[,ee[1,]<0.05]
gg=ee[,ee[1,]<0.01]
cox2=ff[-1,]
write.csv(cox2,'brca_cox2.csv')
