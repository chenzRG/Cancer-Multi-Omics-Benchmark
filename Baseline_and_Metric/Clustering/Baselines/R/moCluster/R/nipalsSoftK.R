nipalsSoftK <-
function(x, maxiter, k) {
  
  if (length(k) < length(x))
    k <- rep(k, length.out = length(x))
  
  t <- svd(do.call("cbind", x))$u[, 1]
  
  regproj <- function(xb, t, k) { 
    pb <- t(xb) %*% t / c(t(t) %*% t)
    pb <- pb/sqrt(sum(pb^2))
    pb <- softK(pb, k)  # soft-thresholding
    tb <- xb %*% pb
    list(tb=tb, pb=pb) # t-score, p-loading
  }
  
  for (i in 1:maxiter) {
    told <- t
    rp <- mapply(SIMPLIFY = FALSE, function(x, y) regproj(x, t, y), x=x, y=k)
    tm <- sapply(rp, "[[", "tb")
    w <- t(tm) %*% t / c(t(t) %*% t)
    w <- w/sqrt(sum(w^2))
    #  w <- w/sum(w)
    t <- tm %*% w
    if (all.equal(c(t), c(told))[1] == TRUE)
      break
    if (i == maxiter)
      cat("  Note: maximum number of iterations was reached, algrithm may not converge.\n")
  }
  
  res <- list(tb = lapply(rp, "[[", "tb"),
              pb = lapply(rp, "[[", "pb"),
              t = t,
              w = w)
  
  return(res)
}
