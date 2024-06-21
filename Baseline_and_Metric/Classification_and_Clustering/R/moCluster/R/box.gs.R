# ==============================================================================
# ==                                                                          ==
# ==                    the box.gs plot                                       ==
# ==                                                                          ==
# ==============================================================================
box.gs <- function(x, gs, moa=NULL, col=1, layout=NULL, plot=TRUE, obs.order=NULL, ...) {
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




