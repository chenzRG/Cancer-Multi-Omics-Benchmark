decompose.gs.ind <- function(x, gs, obs, type=3, nf=2, plot=TRUE,
                   col.data=NULL, col.pc=NULL, legend=TRUE) {
  
  # type = 1 - the data - pc mode
  # type = 2 - the PC -data mode
  # type = 3 - both
  
  if (inherits(x, "mgsa"))
    x <- x@sup
  scl <- lapply(x@score.sep, function (x) x[1:nf])
  m <- sapply(scl, function(x) 
    sapply(x, function(y) mean(y[gs, obs]))) # sum changed here 
  nc <- ncol(m)
  nr <- nrow(m)
  
  if (plot) {
    lgt1 <- lgt2 <- NULL
    if (legend) {
      lgt1 <- names(x@score.pc)[1:nf]
      lgt2 <- names(x@score.data)
    }
    
    if (type == 3)
      layout(matrix(1:2, 1, 2))
    if (type %in%  c(1, 3)) { # data-pc
      ma.d <- max(c(m, colSums(m)))
      mi.d <- min(c(m, colSums(m)))
      barplot(colSums(m), width=nr, space=1/nr, border=NA, ylim=c(mi.d, ma.d), col=col.data)
      barplot(m, beside=TRUE, add=TRUE, border=FALSE, col=col.pc, legend.text=lgt1,
              axes=FALSE)
    } 
    if (type %in% c(2, 3)) { # pc - data
      ma.p <- max(c(m, rowSums(m)))
      mi.p <- min(c(m, rowSums(m)))
      barplot(rowSums(m), width=nc, space=1/nc, border=NA, ylim=c(mi.p, ma.p), col=col.pc)
      barplot(t(m), beside=TRUE, add=TRUE, border=FALSE, col=col.data, legend.text=lgt2,
              axes=FALSE)
    } 
  } 
  return(invisible(m))
}
# ==============================================================================
# ==                                                                          ==
# ==                    the box.gs plot                                       ==
# ==                                                                          ==
# ==============================================================================
box.gs.feature <- function(x, gs, moa=NULL, col=1, layout=NULL, plot=TRUE, obs.order=NULL, ...) {
  # x - either mgsa or moa.sup
  if (inherits(x, "mgsa")) {
    if (!is.null(moa))
      warning("x is an object of class mgsa, moa argument is ignored!")
    data <- x@moa@data
    x <- x@sup
  } else if (inherits(x, "moa.sup")) {
    if (!inherits(moa, "moa"))
      stop("x is an object of class moa.sup, moa argument is required!")
    data <- moa@data
  }

  if (is.null(obs.order))
    obs.order <- 1:ncol(data[[1]])
  col <- rep(col, length(obs.order))[obs.order]

  gsidx <- lapply(x@sup, function(x) as.logical(x[, gs]))
  mats <- mapply(SIMPLIFY=FALSE, function(x, i) x[i, obs.order], x=data, i=gsidx)
# sapply(mats, dim)
  n <- length(mats)
  if (is.null(layout))
    layout <- matrix(1:n, 1, n)
  if (plot) {
    layout(layout)
    for (i in 1:n) {
      if(length(mats[[i]]) == 0)
        plot(0, pch=NA) else
          boxplot(mats[[i]], col=col , ...)
    }
  } else {
    return (mats)
  }
}
