normvec <- function(x) {
  # normalize a vector to unit length
  if (NROW(x) > 1 & NCOL(x) > 1)
    stop("x should be a vector or matrix with 1 row/column.")
  
  if (NCOL(x) > 1)
    prd <- tcrossprod(x) else
      prd <- crossprod(x)
    
  length <- sqrt(c(prd))
  v <- x/length
  attr(v, "length") <- length
  v
}
