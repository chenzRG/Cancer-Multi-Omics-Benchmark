wsvd <- function(X, D1=diag(1, nrow(X)), D2=diag(1, ncol(X))) {
  
  if (is.vector(D1))
    D1 <- diag(D1)
  if (is.vector(D2))
    D2 <- diag(D2)
  i1 <- identical(D1, t(D1))
  i2 <- identical(D2, t(D2))
  if (!(i1 & i2))
    warning("non-symetric distance matrix")
  
  r <- list(D1=D1, D2=D2)
  X <- matpower(D1, 1/2) %*% X
  X <- X %*% matpower(D2, 1/2)
  result <- svd(X)
  r$d <- result$d
  r$u <- matpower(D1, -1/2) %*% result$u
  r$v <- matpower(D2, -1/2) %*% result$v
  
  if (all(r$u[, 1] < 0)) {
    r$u <- r$u * (-1)
    r$v <- r$v * (-1)
  }
  return(r[c("d", "u", "v", "D1", "D2")])
}
