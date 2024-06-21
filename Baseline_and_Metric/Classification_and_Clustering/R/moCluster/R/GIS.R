GIS <- function(x, geneSet, nf=NA, barcol=NA, topN=NA, plot=TRUE, Fvalue=FALSE, ff=NA, cor=FALSE) {
  
  if (!inherits(x, "mgsa"))
    stop("x should be an object of class mgsa!")
  
  fvalue <- function(a, f) {
    n <- length(f)
    ng <- length(unique(f))
    ma <- mean(a)
    ssb <- sum(tapply(a, f, function(x) length(x)*(mean(x)-ma)^2))
    ssw <- sum(tapply(a, f, function(x) sum((x-mean(x))^2)))
    ssb*(n-ng)/ssw/(ng-1)
  }
  
  exprs <- x@moa@data
  rn <- lapply(exprs, rownames)
  rn <- unlist(rn)
  
  if (is.na(nf))
    nf <- length(x@sup@score.pc)
  if (is.na(barcol[1]))
    barcol <- 1:length(exprs)
  
  scor <- x@sup@score[geneSet, , drop=FALSE]
  gsi <- lapply(x@sup@sup, function(x) x[, geneSet])
  coln <- sapply(gsi, function(x) sum(x!=0))
  col_code <- rep(barcol, coln)
  gsi <- unlist(gsi)
  genes_idx <- which(gsi != 0)
  
  if (cor) {
    if (Fvalue)
      warning("cor = TRUE, Fvalue is ignored.")
    xmat <- as.matrix(x@moa@fac.scr[, 1:nf]) %*% t(as.matrix(x@moa@loading)[, 1:nf])
    ef <- cor(t(scor), xmat[, genes_idx], )
    cn <- colnames(ef)
    ef <- c(ef)
    names(ef) <- cn
  } else {
    if (!Fvalue) {
      sdscor <- sd(scor)
    } else if (Fvalue & length(ff)==length(scor)) {
      sdscor <- fvalue(scor, as.factor(ff))
    } else 
      stop("if Fvalue is TRUE, ff need to be defined!")
    gsimat <- sapply(genes_idx, function(i, gsi) {
      gsi[i] <- 0
      return(gsi)}, gsi=gsi)
    coor <- t(gsimat) %*% as.matrix(x@moa@loading)[, 1:nf]
    scor_rm <- as.matrix(x@moa@fac.scr[, 1:nf]) %*% t(coor)
    
    if (!Fvalue) {
      sdscor_rm <- apply(scor_rm, 2, sd)
    } else if (Fvalue & length(ff)==length(scor)) {
      sdscor_rm <- rowFtests(t(scor_rm), fac=as.factor(ff))[, "statistic"]
    }
    names(sdscor_rm) <- rn[genes_idx]
    ef <- -log2(sdscor_rm/sdscor)
  }  
  
  ef <- ef/max(ef)
  oef <- order(ef, decreasing=FALSE)
  ef <- ef[oef]
  ef <- rev(ef)
  
  if (plot) {
    barplot(rev(ef), horiz=TRUE, col=col_code[oef], border=col_code[oef], names.arg = NA)
    legend(x="bottomright", col=barcol, legend=names(coln), pch=15)
  }
  
  dn <- names(x@moa@data)
  fd <- as.factor(rev(col_code[oef]))
  if (!is.null(dn)) levels(fd) <- dn
  r <- data.frame(feature = names(ef), GIS = ef, data=fd)
  r$feature <- as.character(r$feature)
  
  if (!is.na(topN))
    r <- r[1:min(topN, nrow(r)), ]
  return(r)
}


annotate.gs <- function(mgsa, gs) {
  if (length(gs) > 1)
    stop("gs has to be an one element character string or integer")
  rn <- sapply(mgsa@moa@data, rownames)
  r <- lapply(mgsa@sup@sup, function(x) x[, gs] != 0)
  gs <- mapply(SIMPLIFY=FALSE, function(x, g) x[g], x=rn, g=r)
  var <- unique(unlist(gs))
  data <- sapply(gs, function(x) var %in% x) 
  stat <- rowSums(data)
  df <- data.frame(var=var, data=data, stat=stat)
  df <- df[order(df$stat, decreasing=TRUE), ]
  return(df)
}

