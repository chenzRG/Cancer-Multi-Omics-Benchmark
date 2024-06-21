
plotGS <- function(x, axes=1:2, center.only=FALSE, topN=1, data.pch=20, data.col=1, highlight.col = 2,
                          label=NULL, label.cex=1, layout=NULL, ...) {
  
  if (inherits(x, "mgsa")) {
    data <- x@sup
  } else if (inherits(x, "moa.sup")) {
    data <- x
  } else 
    stop("x should be either mgsa or mog.sup class.")

  if (center.only) {
      r <- .plotGS.comb(data, axes=axes, label=label, label.cex=label.cex, 
        topN=topN, data.col=data.col, data.pch=data.pch, highlight.col = highlight.col, ...)
    } else {
      r <- .plotGS.sep(data, axes=axes, label=label, label.cex=label.cex, 
        topN=topN, layout=layout, data.col=data.col, data.pch=data.pch, highlight.col = highlight.col, ...) 
    }
  return(invisible(r))
}

# data is moa.sup object
.plotGS.sep <- function(data, axes, label, label.cex, topN, layout, data.col, data.pch, highlight.col, ...) {

      Nd <- length(data@coord.sep)
      if (length(data.pch)==1) data.pch <- rep(data.pch, Nd)
      if (length(data.col)==1) data.col <- rep(data.col, Nd)
      if (is.null(layout)) lo <- matrix(1:Nd, 1, Nd) else 
        lo <- layout
      layout(lo)

      sg <- lapply(data@coord.sep, .selectTopN, axes, n=topN)
      label <- unique(c(label, unlist(sg)))

      for (i in 1:Nd) {
        x <- data@coord.sep[[i]]
        plot(x[, axes], pch=NA)
        abline(v=0, h=0)
        col <- rep(data.col[i], nrow(x))
        col[rownames(x) %in% label] <- highlight.col
        points(x[, axes], pch=data.pch, col=col, ...)
        text(x=x[label, axes], labels=row.names(x[label, ]), cex=label.cex)
      }
      return(sg)
}

.plotGS.comb <- function(data, axes, label, label.cex, topN, data.col, data.pch, highlight.col, ...) {
  mat <- data@coord.comb
  sg <- .selectTopN(mat, axes, topN)
  label <- unique(c(label, unique(unlist(sg))))
  plot(mat[, axes], pch=NA)
  abline(v=0, h=0)
  col <- rep(data.col[1], nrow(mat))
  col[ rownames(mat) %in% label ] <- highlight.col
  points(mat[, axes], pch=data.pch, col=col, ...)
  text(x=mat[label, axes], labels=row.names(mat[label, ]), cex=label.cex)
  return(sg)
}

.selectTopN <- function(mat, col, n) {
  # give a matrix, select top positve and negative N elements for specific columns.
  if (n != 0) {
      r <- lapply(col, function(i) {
      v <- mat[, i]
      names(v) <- rownames(mat)
      sv <- sort(v)
      lv <- length(v)
      slow <- names(sv[1:n])
      shigh <- names(sv[(lv-n+1):lv])
      list(negN = slow, posN = shigh)
    })
    r <- unlist(r, recursive = FALSE)
    names(r) <- paste(paste("ax", rep(col, each=2), sep=""), names(r), sep=".")
  }
  else r <- NULL
  r
}