concord0 <- function(x, y, ncomp=2, dmod = 1, center = TRUE, scale = FALSE, option = "uniform", 
                    kx = "all", ky = "all", wx = 1, wy = 1, pos = FALSE) {
  
  y <- t(scale(t(y), center = center, scale = scale))
  x <- processOpt(x, center = center, scale = scale, option = option)
  X <- x
  Y <- y
  if (kx == "all") kx <- Inf
  if (ky == "all") ky <- Inf
  
  nmat <- length(X)
  if (is.null(names(X)))
    names(X) <- paste0("X", 1:nmat)
  
  nr <- sapply(X, nrow)
  nc <- unique(sapply(X, ncol))
  if (length(nc) > 1)
    stop("Number of columns in X need to be the same.")
  
  i.sample <- rep(names(X), each = nc)
  i.feature <- rep(names(X), nr)
  
  Ynorm <- scale(t(Y), center = TRUE, scale = TRUE)
  Xnorm <- lapply(X, function(x) scale(t(x), center = TRUE, scale = TRUE))
  Xcat <- do.call("cbind", Xnorm)
  Ynorm.o <- Ynorm
  Xnorm.o <- Xnorm
  
  ys <- yloading <- gls <- bls <- loading <- c()
  xvar <- var <- c()
  
  
  for (f in 1:ncomp) {
    print(f)
    if (f == 1 || dmod != 1 )
      S <- t(Ynorm) %*% Xcat
    if (f == 1)
      S.o <- S
    
    # decom <- propack.svd(S, neig = 1, opts = list(kmax = 20))
    decom <- softSVD(x = S, nf = 1, kv = kx, ku = ky, wv = wx, wu = wy, pos = pos, maxiter = 1000)
    
    xa <- Xcat %*% decom$v[, 1]
    yb <- Ynorm %*% decom$u[, 1]
    var <- c(var, decom$d[1]^2)
    ys <- cbind(ys, yb)
    gls <- cbind(gls, xa)
    yloading <- cbind(yloading, decom$u[, 1])
    loading <- cbind(loading, decom$v[, 1])
    
    xai <- lapply(names(X), function(x) {
      ii <- i.feature == x
      Xcat[, ii] %*% decom$v[ii, 1]
    })
    
    xai.var <- sapply(xai, crossprod)
    xvar <- cbind(xvar, xai.var)
    bls <- cbind(bls, unlist(xai))
    
    if (dmod == 1) {
      # deflation of S, crossprod matrix, the save with SVD directly
      # the classical concordance approach
      S <- S - tcrossprod(decom$u[, 1]) %*% S
      # or the same 
      # Ynorm <- Ynorm - Ynorm %*% tcrossprod(decom$u[, 1])
    } else if (dmod == 2) {
      # deflaltion X using loading of X, as Lafosse & Hanafi 1997
      # but not possible to incorporate with sparse factor
      Xcat <- Xcat - Xcat %*% tcrossprod(decom$v[, 1])
    } else if (dmod == 3) {
      # defaltion Y using its normed score
      Ynorm <- Ynorm - t(t(Ynorm) %*% tcrossprod(normvec(yb)))
    } else if (dmod == 4) {
      # defaltion X and Y using normed score of Y, not suggested
      Ynorm <- Ynorm - t(t(Ynorm) %*% tcrossprod(normvec(yb)))
      Xcat <- Xcat - t(t(Xcat) %*% tcrossprod(normvec(yb)))
    } else if (dmod == 5) {
      # defaltion X and Y using normed score of Y, not suggested
      Ynorm <- Ynorm - t(t(Ynorm) %*% tcrossprod(normvec(yb)))
      Xcat <- Xcat - t(t(Xcat) %*% tcrossprod(normvec(xa)))
    } else {
      stop("unknown dmod")
    }
    
  }
  
  rownames(loading) <- colnames(Xcat)
  
  list(loading.x = loading, 
       loading.y = yloading,
       score.x = bls,
       score.xcomb = gls,
       score.y = ys, 
       loading.x.index = i.feature, 
       score.x.index = i.sample)
}


