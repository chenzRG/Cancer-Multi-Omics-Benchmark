msvd <-
function(x, svd.sol=svd) {
  nd <- length(x)
  nVar <- sapply(x, ncol)
  if (!is.null(names(x)))
    idx <- rep(names(x), times = nVar) else
      idx <- rep(1:length(x), times = nVar)
  nm <- lapply(x, colnames)
  
  tab <- do.call("cbind", x)
  
  dc <- svd.sol(tab)
  
  res <- list()
  res$t <- dc$u[, 1, drop=FALSE] * dc$d[1]
  pb <- split(dc$v[, 1, drop=FALSE], idx)
  for (i in names(x)) names(pb[[i]]) <- colnames(x[[i]]) 
  res$pb <- lapply(pb, function(x) as.matrix(x/sqrt(sum(x^2))))
  
  res$tb <- mapply(SIMPLIFY = FALSE, function(m, v) {
    m %*% v
  }, m=x, v=res$pb[names(x)])
  
  tm <- do.call("cbind", res$tb)
  res$w <- t(tm) %*% res$t / c(t(res$t) %*% res$t)
  rownames(res$w) <- names(x)
  
  res <- res[c("tb", "pb", "t", "w")]
  res$tb <- res$tb[names(x)]
  res$pb <- res$pb[names(x)]
  return(res)
}
