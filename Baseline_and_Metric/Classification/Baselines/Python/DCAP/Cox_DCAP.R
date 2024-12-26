library("xgboost")
library("ipred")
library("survival")
library("survivalROC")
library("glmnet")
library('kernlab')
library('plyr')
library('caret')

setwd('D:/brca')

kx2=read.csv('data_cox2.csv',row.names= 1)
sdata=read.csv('brca_go.csv',row.names= 1)
data1=t(sdata)

j=6
set.seed(j)
k=10
folds <- createFolds(as.data.frame(t(kx2)),k)
results<-matrix(nrow = k, ncol=6)
im=matrix(nrow=500,ncol=30)
eva<-matrix(nrow = nrow(kx2), ncol=k)
count=0
for(i in 1:k){
  testset <- kx2[folds[[i]],] 
  trainset <- kx2[-folds[[i]],] 
  te_xg=data1[folds[[i]],] 
  tr_xg=data1[-folds[[i]],] 
  x=trainset[,-c(1:2)]
  x=as.matrix(x)
  tc=trainset$time
  tc[tc==0]=0.001
  y=Surv(tc,trainset$status)
  cv.fit<-cv.glmnet(x,y,family="cox",maxit=10000,alpha=0,nfold=5)
  fit<-glmnet(x,y,family="cox",alpha=0)
  tt=predict(fit,x,s=cv.fit$lambda.min)
  x2=testset[,-c(1:2)]
  x2=as.matrix(x2)
  tc2=testset$time
  tc2[tc2==0]=0.001
  y2=Surv(tc2,testset$status)
  tt2=predict(fit,x2,s=cv.fit$lambda.min)
  y_xgtr=tt;
  x_xgtr=tr_xg
  x_xgtr=as.matrix(x_xgtr)
  x_xgte=te_xg
  x_xgte=as.matrix(x_xgte)
  bst <- xgboost(x_xgtr, y_xgtr, 
                 max_depth =3, eta =0.25, nrounds =6, min_child_weight=2,
                 objective = "reg:linear",eval_metric = "rmse")
  pred <- predict(bst, x_xgte)
  eva[(count+1):(count+nrow(tt2)),1]<-tc2
  eva[(count+1):(count+nrow(tt2)),2]<-testset$status
  eva[(count+1):(count+nrow(tt2)),3]=tt2
  eva[(count+1):(count+nrow(tt2)),4]=pred
  count=count+nrow(tt2)
  
}

tc3=eva[,1]
st3=eva[,2]
rk3=eva[,3]
rk4=eva[,4]

y3=Surv(tc3,st3)
ci_Cox=survConcordance(formula = y3~ rk3)
ci_Cox$concordance
ci_XGB=survConcordance(formula = y3~ rk4)
ci_XGB$concordance


