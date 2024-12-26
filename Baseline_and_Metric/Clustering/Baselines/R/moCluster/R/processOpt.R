processOpt <-
function(x, center=TRUE, scale=FALSE, option = c("lambda1", "inertia", "uniform")) {
  
  opt <- match.arg(option)  
  
  if (is.null(names(x)))
    names(x) <- paste("data", 1:length(x), sep = "_")

  x <- lapply(x, scale, center, scale)
  if (opt == "lambda1") {
    w <- sapply(x, function(xx) 1/svd(xx)$d[1])
  } else if (opt == "inertia") {
    w <- sapply(x, function(xx) 1/sqrt(sum(xx^2)))
  } else if (opt == "uniform") {
    w <- rep(1, length(x))
  }
  mapply(SIMPLIFY = FALSE, function(xx, ww) xx*ww, xx=x, ww=w)
}
