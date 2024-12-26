
.ks.pval <- function(moa, sup, ks.B, A, factors, mc.cores = getOption("mc.cores", 2L)) {
  
  rmat <- as.matrix(moa@loading[, factors]) %*% t(as.matrix(moa@fac.scr[, factors])) * A
  sup <- do.call("rbind", sup)
  supIdx <- apply(sup, 2, function(x) which(x != 0))
  if (is.matrix(supIdx))
    supIdx <- as.data.frame(supIdx)
  v <- .gsva.ssgsea(rmat, geneSets = supIdx)
  
  ll <- mclapply(1:ks.B, function(x) {
    if (x%%50==0)
      cat(paste("iter", x, "is finished.\n"))
    mm <- apply(rmat, 2, sample, replace=TRUE)
    c(.gsva.ssgsea(mm, geneSets = supIdx))
  }, mc.cores = mc.cores)
  cp <- do.call("cbind", ll)
  pv <- sapply(1:length(v), function(i) {
    p1 <- 2*(1 + sum(cp[i, ] > v[i]))/(ks.B+1)
    p2 <- 2*(1 + sum(cp[i, ] < v[i]))/(ks.B+1)
    pmin(p1, p2)
  })
  pv[pv > 1] <- 1
  pv <- matrix(pv, nrow = nrow(v), 
               dimnames = list(rownames(v), colnames(v)))
  pv
}

.gsva.rndWalk <- function(gSetIdx, geneRanking, j, R, alpha) {
  # this function is copyied from GSVA package
  indicatorFunInsideGeneSet <- match(geneRanking, gSetIdx)
  indicatorFunInsideGeneSet[!is.na(indicatorFunInsideGeneSet)] <- 1
  indicatorFunInsideGeneSet[is.na(indicatorFunInsideGeneSet)] <- 0
  stepCDFinGeneSet <- cumsum((abs(R[geneRanking, j]) * 
                                indicatorFunInsideGeneSet)^alpha) /
    sum((abs(R[geneRanking, j]) *
           indicatorFunInsideGeneSet)^alpha)
  stepCDFoutGeneSet <- cumsum(!indicatorFunInsideGeneSet) /
    sum(!indicatorFunInsideGeneSet)
  walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet
  sum(walkStat) 
}

.gsva.ssgsea <- function(X, geneSets, alpha=0.25) {
  # this function is modified from GSVA package
  p <- nrow(X)
  n <- ncol(X)
  R <- apply(X, 2, function(x) rank(x, ties.method = "min"))
  es <- sapply(1:n, function(j, R, geneSets, alpha) {
    geneRanking <- order(R[, j], decreasing=TRUE)
    es_sample <- sapply(geneSets, .gsva.rndWalk, geneRanking, j, R, alpha)
    unlist(es_sample)
  }, R, geneSets, alpha)
  if (length(geneSets) == 1)
    es <- matrix(es, nrow=1)
  rownames(es) <- names(geneSets)
  colnames(es) <- colnames(X)
  es
}

