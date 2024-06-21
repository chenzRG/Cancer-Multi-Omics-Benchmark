
concordanceNetworkNMI = function(Wall,C) {
	# Given a list of affinity matrices, Wall, the number of clusters C, return a matrix containing the NMIs
# between cluster assignments made with spectral clustering.


#For example, if the input "wall" contains two networks, the output is a 3x3 matrix of NMIs. Note here NMI is a metric to measure the fitness of two clustering results. if NMI = 1, it indicates two clusters are identical; if NMI = 0, it indicates two #clusters are totally differently, independent. 
# the output looks like this:
#                                 fused network     network1        newtwork2
#fused network                 1                       0.7                      0.55
#network    1                   0.72                      1                         0.67
#network 2        0.55        0.67         1

#The output above means that, the fused network is more similar to network 1 than network 2.  And network 1 and network 2 has a similarity 0.67. 
#sometimes, if the inputed networks have low similarities, which means they are contradicting each other, fusing them may not be a good idea.


  # Calculate the fused network
  # C is the number of clusters.
  LW = length(Wall)
   # Get the cluster labels for each of the networks
  labels = lapply(Wall, function(x) spectralClustering(x, C))
  # Calculate the NMI between each pair clusters
  NMIs = matrix(NA, LW, LW)
  for (i in 1:LW) {
    for (j in 1:LW) {
      NMIs[i, j] = calNMI(labels[[i]], labels[[j]]);
    }
  }
  
  return(NMIs)
}
