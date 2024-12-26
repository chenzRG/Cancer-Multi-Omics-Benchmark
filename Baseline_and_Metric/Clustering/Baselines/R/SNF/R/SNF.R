SNF <- function(Wall,K=20,t=20) {
	
	###This function is the main function of our software. The inputs are as follows:
    # Wall : List of affinity matrices
    # K : number of neighbors
    # t : number of iterations for fusion
    
    ###The output is a unified similarity graph. It contains both complementary information and common structures from all individual network. 
    ###You can do various applications on this graph, such as clustering(subtyping), classification, prediction.
    
    LW = length(Wall)
    normalize <- function(X) X / rowSums(X)
    # makes elements other than largest K zero
    
    
    newW <- vector("list", LW)
    nextW <- vector("list", LW)
    ###First, normalize different networks to avoid scale problems.
    for( i in 1: LW){
      Wall[[i]] = normalize(Wall[[i]]);
      Wall[[i]] = (Wall[[i]]+t(Wall[[i]]))/2;
    }
    
    ### Calculate the local transition matrix.
    for( i in 1: LW){
      newW[[i]] = (.dominateset(Wall[[i]],K))
    }
    
    # perform the diffusion for t iterations
    for (i in 1:t) {
      for(j in 1:LW){
        sumWJ = matrix(0,dim(Wall[[j]])[1], dim(Wall[[j]])[2])
        for(k in 1:LW){
          if(k != j) {
            sumWJ = sumWJ + Wall[[k]]
          }
        }
          nextW[[j]] = newW[[j]] %*% (sumWJ/(LW-1)) %*% t(newW[[j]]);
      }
      ###Normalize each new obtained networks.
      for(j in 1 : LW){
      	      
                Wall[[j]] = nextW[[j]] + diag(nrow(Wall[[j]]));
                Wall[[j]] = (Wall[[j]] + t(Wall[[j]]))/2;
                }
   }
    
    # construct the combined affinity matrix by summing diffused matrices
    W = matrix(0,nrow(Wall[[1]]), ncol(Wall[[1]]))
    for(i in 1:LW){
      W = W + Wall[[i]]
    }
    W = W/LW;
    W = normalize(W);
    # ensure affinity matrix is symmetrical
    W = (W + t(W)+diag(nrow(W))) / 2;
    
    return(W)  
  }
