mbpca <- 
  function (x, ncomp, method, k = "all", center = TRUE, 
    scale = FALSE, option = "uniform", maxiter = 1000, 
    moa = TRUE, verbose = TRUE, svd.solver = c("svd", "fast.svd", "propack"),
    k.obs = "all", w = NA, w.obs = NA,
    unit.p = FALSE, unit.obs = FALSE, pos = FALSE) {
    
    method <- match.arg(method, c("globalScore", "blockScore", "blockLoading"))
    x <- lapply(x, t)

    call <- match.call()
    x <- processOpt(x, center = center, scale = scale, option = option)
    nc <- sapply(x, ncol)
    keepAllb <- k[1] == "all"
    keepAllt <- k.obs[1] == "all"
    
    prddata <- lapply(x, t)
    ssl <- match.arg(svd.solver[1], c("svd", "fast.svd", "propack"))
    svdf <- switch(ssl, 
                   "svd" = svd, 
                   "fast.svd" = fast.svd, 
                   "propack" = function(X) propack.svd(X, neig = 1, opts = list(kmax = 20)))
    for (i in 1:ncomp) {
      if (verbose) 
        cat(paste("calculating component ", i, " ...\n", sep = ""))
      if (keepAllb & keepAllt) 
        r <- msvd(x, svd.sol = svdf)
      else {
        
        if (is.na(w)[1])
          w <- 1
        if (is.na(w.obs)[1])
          w.obs <- 1
        
        if (keepAllb == "all")
          keepAllb <- Inf
        if (keepAllt == "all")
          keepAllt <- Inf
        if (length(k) < length(x))
          k <- rep(k, length.out = length(x))
        if (length(k.obs) < length(x))
          k.obs <- rep(k.obs, length.out = length(x))
        r <- biSoftK(x, maxiter = maxiter, kp = k, kt = k.obs, unit.pb = unit.p, 
          unit.tb = unit.obs, weight.p = w, weight.t = w.obs, pos = pos)
      }
      x <- deflat(x, r$t, r$tb, r$pb, method)
      if (i == 1) {
        res <- r 
      } else {
        res$t <- cbind(res$t, r$t)
        res$w <- cbind(res$w, r$w)
        res$tb <- mapply(cbind, res$tb, r$tb, SIMPLIFY = FALSE)
        res$pb <- mapply(cbind, res$pb, r$pb, SIMPLIFY = FALSE)
      }
    }
    if (moa) 
      res <- toMoa(prddata, res, call = call)
    return(res)
  }
