decompose.gs.group <-
function(x, gs, group, decomp = "data", nf=2, x.legend="bottomleft", y.legend=NULL, plot=TRUE, 
  main=NULL, ...) {
  
  # function only used inside function
  sedata <- function(ob) {
    m <- sapply(scl, function(x) 
      rowSums(sapply(x, function(y) y[gs, ob])))
    se <- rowSds(t(m))/sqrt(nrow(m))
    return(se)
  }
  
  sepc <- function(ob) {
    m <- sapply(scl, function(x) 
      sapply(x, function(y) y[gs, ob]))
    l <- split(m, rep(1:nf, rep(length(ob), nf)))
    l <- lapply(l, function(x) x[x!= 0])
    se <- sapply(l, sd)/sqrt(length(l[[1]]))
    names(se) <- names(scl[[1]])
    return(se)
  }
  # done
  
  if (inherits(x, "mgsa"))
    x <- x@sup
  scl <- lapply(x@score.sep, function (x) x[1:nf])
  
  if (is.list(group))
    cls <- group else {
      if (is.null(names(group)))
        names(group) <- paste("V", 1:length(group), sep = "")
        cls <- split(1:length(group), group)
    }
    
  gsm <- lapply(cls, function(ob) 
    decompose.gs.ind(x = x, plot = FALSE, gs=gs, obs = ob, nf = nf))
  
  if (decomp == "data") {
   s <- sapply(gsm, colSums) # by data
   sed <- sapply(cls, sedata)
  } else if (decomp == "pc") {
    s <- sapply(gsm, rowSums) # by pc
    sed <- sapply(cls, sepc)
  } else
    stop("unknown setting of decomp, decomp should be either data or pc!")
  
  if (plot) {
    ci.l <- s-1.96*sed
    ci.u <- s+1.96*sed
    col <- gray.colors(nrow(s))
    u <- barplot2(s, beside = TRUE, plot.ci = TRUE, ci.l = ci.l, ci.u = ci.u, col=col, 
      ylab = paste(decomp, "-wise decomposed gene set scores", sep=""), main=main)
    legend(x = x.legend, y = y.legend, legend = rownames(s), col=col, pch=15)
  }
  return(invisible(list(decomp.mean=s, decomp.se=sed)))
}
