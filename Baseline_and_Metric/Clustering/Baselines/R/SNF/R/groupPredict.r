groupPredict <- function(train,test,groups,K=20,alpha=0.5,t=20,method=1){

###This function is used to predict the subtype of new patients.
#train and test have the same number of view and the same number of columns
# group is the label for the train data
# K, alpha, t are the prameters for SNF. 
#K is the number of neighbors
#alpha is the hyperparameter used in constructing similarity network
# t is the number of iterations
#method is a indicator of which method to use to predict the label. method = 0 means to use local and global consistency; method = 1 means to use label propagation.
Wi= vector("list", length=length(train));

for (i in 1:length(train)){
view= standardNormalization(rbind(train[[i]],test[[i]]));
Dist1 = dist2(view, view);
Wi[[i]] = affinityMatrix(Dist1, K, alpha);
}

W = SNF(Wi,K,t);
Y0=matrix(0,nrow(view), max(groups));
for (i in 1:length(groups)) Y0[i,groups[i]]=1;
Y=.csPrediction(W,Y0,method);
newgroups=rep(0,nrow(view));
for (i in 1:nrow(Y)) newgroups[i]=which(Y[i,]==max(Y[i,]));

return (newgroups);
}
