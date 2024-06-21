
moaScore <- function(moa) moa@fac.scr

moaCoef <- function(moa) {
  mm <- moa@loading
  m <- as.data.frame(mm)
  nd <- length(moa@data)
  fac <- rep(1:nd, moa@tab.dim[1, ])
  ms <- split(m, fac)
  
  r <- lapply(ms, function(x) {
    var <- rownames(x)
    lapply(x, function(y) {
      names(y) <- var
      neg <- sort(y[y<0], decreasing = FALSE)
      pos <- sort(y[y>0], decreasing = TRUE)
      
      neg <- data.frame(id=names(neg), coef=neg)
      pos <- data.frame(id=names(pos), coef=pos)
      list(neg=neg, pos=pos)
      })
  })
  names(r) <- names(moa@data)
  list(coefMat = mm, nonZeroCoef = unlist(lapply(r, unlist, recursive=FALSE), recursive = FALSE))
}

