distMoa <-
function(x, nf=NA, tol=1e-5, method = "euclidean", diag = FALSE, upper = FALSE, p = 2) {
  if (is.na(nf))
    nf <- Inf
  if (inherits(x, "moa")) {
    nfi <- x@eig > tol
    x <- moaScore(x)
  }
  if (nf > ncol(x) | nf > sum(nfi)) {
    nf <- min(ncol(x), sum(nfi))
    cat(paste("nf set to ", nf, ".\n", sep = ""))
  }
  nfi[-(1:nf)] <- FALSE
  x <- x[, nfi]
  dist(x, method = method, diag = diag, upper = upper, p = p)
}
