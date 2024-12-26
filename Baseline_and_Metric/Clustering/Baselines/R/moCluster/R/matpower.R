matpower <- function(x, n, nf=min(dim(x)), tol=1e-7) {
  
  if (inherits(x, "matrix")) {
    s <- svd(x, nu = nf, nv = nf)
    m <- sum(s$d > tol)
    nf <- min(nf, m)
    x <- s$u[, 1:nf] %*% diag(s$d[1:nf]^n) %*% t(s$v[, 1:nf])
  } else 
    stop("x need to be a matrix object!")
  return(x)
}

