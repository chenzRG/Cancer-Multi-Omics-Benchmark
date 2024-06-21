dist2 <- function(X,C) {
  ndata = nrow(X)
  ncentres = nrow(C)

  sumsqX = rowSums(X^2)
  sumsqC = rowSums(C^2)
  
  XC = 2 * (X %*% t(C))
  
  res = matrix(rep(sumsqX,times=ncentres),ndata,ncentres) + t(matrix(rep(sumsqC,times=ndata),ncentres,ndata)) - XC
  res[res < 0] = 0
  return(res)
}
