sup.moa <- function(X, sup, nf = 2, factors = NULL, 
  ks.stat=FALSE, ks.B = 1000, ks.cores = NULL, p.adjust.method = "fdr") {

  if (is.null(nf) & is.null(factors))
    stop("nf or factors need to be specified.")

  if (!is.null(factors)) {
    if (!is.null(nf))
      message("factors is given, nf will be ignored.")
    pcomp <- factors
  } else if (!is.null(nf)) 
    pcomp <- 1:nf

  # sup, nf
  N <- length(sup)
  fs <- X@fac.scr 
  
  # sup data and moa data have the same order
  repr <- sapply(X@data, dim)[1, ]
  nn <- names(X@data)
  load <- split(X@loading, f=rep(nn, repr))
  load <- load[nn]
  
  cols <- sapply(sup, colSums)
  if (!is.matrix(cols))
    cols <- matrix(cols, nrow = length(sup))
  w <- rowSums(cols)

  if (any(w == 0)) {
    message("unrelated gene sets are detected and will be removed")
    sup <- lapply(sup, function(x) x[, w > 0])
  }
    
  # normsup <- lapply(sup, function(x, w) {
  #   sweep(x, 2, w, "/")
  # }, w=w)
  normsup <- sup # change here

  GSCoordinate_sep <- mapply(SIMPLIFY=FALSE, function(load, sup, A) {
    a <- t(sup * A) %*% as.matrix(load[, pcomp, drop=FALSE])
    colnames(a) <- paste("PC", pcomp, sep="")
    return(a)
  }, load=load, sup=normsup, A = split(X@w.data, names(X@w.data))[nn])

  GSCoordinate_comb <- Reduce("+", GSCoordinate_sep)
  
  contribution <- lapply(GSCoordinate_sep, function(supcor, score) {
    a <- lapply(1:length(pcomp), function(i) {
      r <- outer(supcor[, i], score[, pcomp[i]])
      colnames(r) <- rownames(score)
      return(r)
    })
    a[is.na(a)] <- 0
    names(a) <- paste("PC", pcomp, sep="")
    return(a)
  }, score=fs)
  
  contribution_dataset <- lapply(contribution, function(x) {
    Reduce("+", x)
  })
  
  contribution_pc <- lapply(1:length(pcomp), function(i, cont) {
    a <- lapply(cont, function(x) x[[i]])
    Reduce("+", a)
  }, cont=contribution)
  names(contribution_pc) <- paste("PC", pcomp, sep="")
  
  contribution_total <- Reduce("+", contribution_dataset) 

  csup <- do.call("rbind", sup)
  if (!ks.stat) {
    pmat <- .signifGS(X=X, sup=csup, A = X@w.data, 
      score=contribution_total, factors=pcomp)
    attr(pmat, "method") <- "zscore"
    } else {
      if (is.null(ks.cores)) 
        ks.cores <- getOption("mc.cores", 2L) 
      cat("running bootstrapping for p values of KS.stat ...\n")
      pmat <- .ks.pval(X, sup, ks.B=ks.B, A = X@w.data, factors=pcomp, mc.cores = ks.cores)
      attr(pmat, "method") <- "KS.stat"
    }
  
  pmatadj <- matrix(p.adjust(pmat, method = p.adjust.method), nrow(pmat), ncol(pmat))
  attr(pmatadj, "method") <- p.adjust.method

  res <- new("moa.sup", 
    sup = sup,
    coord.comb = GSCoordinate_comb,
    coord.sep = GSCoordinate_sep,
    score = contribution_total,
    score.data = contribution_dataset,
    score.pc = contribution_pc,
    score.sep = contribution,
    p.val = pmat,
    p.val.corrected = pmatadj
    )
  return(res)
}


.signifGS <- function(X, sup, A, score, factors) {
  
  # define function 
  ff <- function(x, n, score, infinite=FALSE) {
    lx <- length(x) 
    if (infinite)
      sf <- 1 else
        sf <- sqrt((lx-n)/(lx-1))
    sum_sd <- sf * sd(x)/sqrt(n) * n
    sum_mean <- mean(x) * n
    pp <- abs(pnorm(score, mean = sum_mean, sd = sum_sd))
    2 * rowMin(cbind(pp, 1-pp))
  }
  # reconstuct matrix using nf PCs
  U <- as.matrix(X@loading[, factors, drop=FALSE]) 
  D <- diag(sqrt(X@eig[factors]), nrow = length(factors))
  V <- as.matrix(X@eig.vec[, factors, drop=FALSE])
  rec <- (U %*% D %*% t(V)) * A
  # the number of feature in each GS
  supn <- colSums(sup != 0)
  # calculate the P value
  pmat <- sapply(1:ncol(rec), function(i) ff(rec[, i], supn, score = score[, i]))
  if (!is.matrix(pmat))
    pmat <- matrix(pmat, ncol = ncol(rec))
  
  colnames(pmat) <- colnames(score)
  rownames(pmat) <- rownames(score)
  return(pmat)
}