pairwise.rv <- function(data.list, match="col") {
	if (match %in% "row")
		data.list <- lapply(data.list, t)
  ms <- sapply(data.list, function(x) {
    x <- c(crossprod(as.matrix(x)))
    x <- x/sqrt(sum(x^2))})
  m <- crossprod(ms)
  colnames(m) <- rownames(m) <- names(data.list)
  return(m)
}