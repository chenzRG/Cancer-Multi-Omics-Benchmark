estimateNumberOfClustersGivenGraph <- function(W, NUMC=2:5) {
  
  #   This function estimates the number of clusters given the two huristics
  #   given in the supplementary materials of our nature method paper
  #   W is the similarity graph
  #   NUMC is a vector which contains the possible choices of number of
  #   clusters.
  #   
  #   
  #   K1 is the estimated best number of clusters according to eigen-gaps
  #   K12 is the estimated SECOND best number of clusters according to eigen-gaps
  #   
  #   K2 is the estimated number of clusters according to rotation cost
  #   K22 is the estimated SECOND number of clusters according to rotation cost
  #   
  #   an example would be [K1, K2, K12,K22] = Estimate_Number_of_Clusters_given_graph(W,
  #                                                                                   [2:5]);
  #   
  #   Note that this function can only give an estimate of the number of
  #   clusters. How to determine the "OPTIMAL" number of clusters, is still an
  #   open question so far. 
  
  if (min(NUMC) == 1) {
    warning('Note that we always assume there are more than one cluster.');
    NUMC = NUMC[NUMC > 1]  
  }
  
  W = (W + t(W))/2
  diag(W) = 0
  
  if (length(NUMC) > 0) {
    degs = rowSums(W)

    
    # compute unnormalized Laplacian
    
    degs[degs == 0] = .Machine$double.eps    
    D = diag(degs)    
    L = D - W
    Di = diag(1 / sqrt(degs))
    L = Di %*% L %*% Di
    
    # compute the eigenvectors corresponding to the k smallest
    # eigs$valuess
    eigs = eigen(L)
    eigs_order = sort(eigs$values, index.return=T)$ix
    eigs$values = eigs$values[eigs_order]
    eigs$vectors = eigs$vectors[, eigs_order]
    eigengap = abs(diff(eigs$values))
    eigengap = eigengap * (1 - eigs$values[1:length(eigs$values) - 1] ) / (1 - eigs$values[2:length(eigs$values)])

    quality = list()
    for (c_index in 1:length(NUMC)) {
      ck = NUMC[c_index]
      UU = eigs$vectors[, 1:ck]
      EigenvectorsDiscrete <- .discretisation(UU)[[1]]
      EigenVectors = EigenvectorsDiscrete^2
      
      # MATLAB: sort(EigenVectors,2, 'descend');
      temp1 <- EigenVectors[do.call(order, lapply(1:ncol(EigenVectors), function(i) EigenVectors[, i])), ]
      temp1 <- t(apply(temp1, 1, sort, TRUE))  
      
      quality[[c_index]] = (1 - eigs$values[ck + 1]) / (1 - eigs$values[ck]) * 
        sum( sum( diag(1 / (temp1[, 1] + .Machine$double.eps) ) %*% temp1[, 1:max(2, ck-1)] ))
    }
    
    t1 <- sort(eigengap[NUMC], decreasing=TRUE, index.return=T)$ix
    K1 = NUMC[t1[1]]
    K12 = NUMC[t1[2]]
    t2 <- sort(unlist(quality), index.return=TRUE)$ix
    K2 <- NUMC[t2[1]]
    K22 <- NUMC[t2[2]]    
  }
  
  return (list(K1, K12, K2, K22))
}
