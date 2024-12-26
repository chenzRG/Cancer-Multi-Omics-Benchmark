# plot.moa

# value options:
#   eig - plot the eigen values
#     ... could be:
#       type=1 - the type of plot
#       axes=NULL - the axes selected to plot
#       n=NULL - n eigenvalues to be drawn
#       tol=1e-5 - the tolerance of eigenvalue, eigenvalues lower than this value wont be considered.
#       legend=NULL - legend to put
#       col=NULL - the color of each partial eigenvalue
#       lty=1 - the line type used in the matplot, when type =4, used
#       pch=NULL - the pch to draw 2D partial eigen plot, when type = 5 used
#       lg.x="topright" - the position of legend
#       lg.y=NULL - poistion argument passed to legend(...)
#       ... - other arguemnts passed to functions, see below

#     for:
#       type 1: the eigen value
#         ... are passed to barplot
#       type 2: barplot show, partial eigenvalue, beside=FALSE
#         ... are passed to barplot
#       type 3: barplot show, partial eigenvalue, beside =TRUE
#         ... are passed to barplot
#       type 4: matplot show
#         ... are passed to matplot
#       type 5: the two dimensional plot, axes need to be specified
#         ... are passed to heatmap
  
#   tau - the same with eig, but in the percentage view
#     ... could be (same with eig, but in the percentage)

#   obs - the observation
#     ... could be:
#       axes=1:2 - which axes should be draw
#       type=1 - which type, see below
#       data.pch=20 - the pch of dataset, if type=1, the first one is used
#       col=1 - the color of observations, recycled used by data.frame
#       label=FALSE - should be labeled?
#       lg.x="topright" - position of legend
#       lg.y=NULL - position of legend
#       xlim=NULL - the xlimit
#       ylim=NULL - the ylimit
#       label.cex=1 - the cex of text
#       ...
    
#     for:
#       type 1: the center points draw
#         ... passed to points
#       type 2: the separate factor scores linked by lines
#         ... passed to points
#   var - the separate gene view, layout can be specified
#   RV - the heatmap of RV coefficient

setMethod("plot", signature=c("moa", "missing"), function(x, value, type=1,
  axes=NULL, n=NULL, tol=1e-5, legend=NULL, col=NULL, lty=1, 
  pch=NULL, lg.x="topright", lg.y=NULL, xlim=NULL, ylim=NULL, 
  data.pch=20, label=FALSE, label.cex=1, layout=NULL, ...) {

  if (value %in% c("eig", "eigenvalue")) {
    .plot.eig(x, type=type, axes=axes, n=n, tol=tol, 
      legend=legend, col=col, lty=lty, pch=pch, lg.x=lg.x, lg.y=lg.y, ...)
  } else if (value %in% c("tau", "Perc.eig")) {
    .plot.tau(x, type=type, axes=axes, n=n, tol=tol, 
      legend=legend, col=col, lty=lty, pch=pch, lg.x=lg.x, lg.y=lg.y, ...)
  } else if (value %in% c("obs", "observation")) {
    if (is.null(col))
      col <- "gray25"
    if (is.null(axes))
      axes <- 1:2
    .plot.obs(x, axes=axes, type=type, data.pch=data.pch, col=col,
      label=label, lg.x="topright", lg.y=lg.y, xlim=xlim, ylim=ylim, label.cex=label.cex, ...)
  } else if (value %in% c("var", "variable", "feature")) {
    if (is.null(col))
      col <- "gray25"
    if (is.null(axes))
      axes <- 1:2
    .plot.var(x, axes=axes, layout=layout, col=col, ...)
  } else if (value %in% c("RV", "rv")) {
    if (is.null(col))
      col <- heat.colors(12)
    .heatmap.rv(x, ...)
  } else
    stop ("unknow value selected.")
  
})


# ==============================================================================
# ==                                                                          ==
# ==                     plot.obs                                             ==
# ==                                                                          ==
# ==============================================================================

# type = 1, only plot the center
# type = 2, separate data link by lines

.plot.obs <- function(moa, axes=1:2, type=1, data.pch=20, col=1,
  label=FALSE, lg.x="topright", lg.y=NULL, xlim=NULL, ylim=NULL, label.cex=1, ...) {
  
  if (type == 1) {
    if (length(data.pch) > 1)
      warning("for type 1 plot, only the first data.pch will be used.")
    plot(moa@fac.scr[, axes], pch=NA, xlim=xlim, ylim=ylim)
    abline(v=0, h=0)
    if (length(data.pch) > 1)
      warning("Type 1 plot only do not distinguish datasets, so the 
        first elements in data.pch is used")
    points(moa@fac.scr[, axes], pch=data.pch[1], col=col, ...)

  } else if (type == 2) {

    ndata <- length(moa@partial.fs)
    didx <- sapply(moa@partial.fs, nrow)
    
    if (length(data.pch) == 1)
      data.pch=rep(data.pch, ndata)
    
    cpartfi <- do.call("rbind", moa@partial.fs)
    plot(cpartfi[, axes], pch=NA, xlim=xlim, ylim=ylim)  
    legend(x=lg.x, y=lg.x, legend=names(moa@partial.fs), pch=data.pch)
    abline(v=0, h=0)
    center <- moa@fac.scr[, axes]/ndata 
    for (x in moa@partial.fs)
      segments(x[, axes[1]], x[, axes[2]], center[, 1], center[, 2], col=col)  
    points(cpartfi[, axes], pch=rep(data.pch, didx), 
      col=rep(col, ndata), ...)
  } else {
    stop("unknown type selected")
  }

  if (label) 
    text(moa@fac.scr[, axes[1]], moa@fac.scr[, axes[2]], rownames(moa@fac.scr), 
      cex=label.cex)
}

# ==============================================================================
# ==                                                                          ==
# ==                     plot.variables/features                              ==
# ==                                                                          ==
# ==============================================================================

.plot.var <- function(moa, axes=1:2, layout=NULL, col="gray25", ...) { 

  tabidx <- sapply(moa@data, dim)[1, ]
  loadings <- split(moa@loading, rep(names(tabidx), tabidx))
  loadings <- loadings[names(tabidx)]
  if (is.null(layout))
    layout <- matrix(1:length(tabidx))
  layout(layout)
  t <- lapply(names(loadings), function(x) {
    mmt <- loadings[[x]]
    plot(mmt[, axes], pch=NA, main = x)
    abline(v=0, h=0)
    points(mmt[, axes], pch=20, col=col, ...)
  })
}

# ==============================================================================
# ==                                                                          ==
# ==               contribution of eigen value (tabls)                        ==
# ==                                                                          ==
# ==============================================================================

# for type 1 to 3, ... are argumens passed to barplot
# for type equals 4, ... are arguments passed to matplot
# for type equals 5, ... are arguments passed to plot

.plot.eig <- function(moa, type=1, axes=NULL, n=NULL, tol=1e-7, 
                      legend=NULL, col=NULL, lty=1, pch=NULL, lg.x="topright", lg.y=NULL, ...) {
  
  if (length(type) != 1)
    warning("only the first value in type is used!")
  if (!type[1] %in% 1:5)
    stop("Unknown type sepecified!")
  
  s <- sum(moa@eig > tol)
  if (!is.null(n))
    s <- min(n, s)  
  ei <- (moa@eig/moa@eig[1])[1:s]
  m <- as.matrix(moa@partial.eig[, 1:s])
  
  if (is.null(legend))
    legend <- rownames(m)
  
  if (type == 1) {
    barplot(moa@eig[1:s], col=col, ...)
    if (!is.null(axes))
      warning("type 1 renders barplot, axes argument are ignored")
  } else if( type == 2) { 
    barplot(m, legend.text=legend, plot=TRUE, col=col, ...)
    if (!is.null(axes))
      warning("type 2 renders barplot, axes argument are ignored")
  } else if (type == 3) {
    barplot(m, beside=TRUE, legend.text=legend, col=col,  ...)
    if (!is.null(axes))
      warning("type 3 renders barplot, axes argument are ignored")
  }  else if (type == 4) {
    if (is.null(col)) col <- 1
    matplot(t(m), type="l", col=col, lty=lty, ...)
    legend(x=lg.x, y=lg.y, legend=legend, col=col, lty=lty)
    if (!is.null(axes))
      warning("type 4 renders lines, axes argument are ignored")
  } else if (type == 5) {
    if (!inherits(axes, "numeric") | length(axes) != 2)
      stop("axes must consists of two integers")
    if (is.null(col)) col <- 1
    plot(m[, axes], col=col, pch=pch, ...)
    legend(x=lg.x, y=lg.y, legend=legend, pch=pch, col=col)
  }
}


# ==============================================================================
# ==                                                                          ==
# ==               contribution of tabs                                       ==
# ==                                                                          ==
# ==============================================================================

# for type 1 to 3, ... are argumens passed to barplot
# for type equals 4, ... are arguments passed to matplot (except type)
# for type equals 5, ... are arguments passed to plot

.plot.tau <- function(moa, type=1, axes=NULL, n=NULL, tol=1e-5, 
  legend=NULL, col=NULL, lty=1, pch=NULL, lg.x="topright", lg.y=NULL, ...) {
  
  if (length(type) != 1)
    warning("only the first value in type is used!")
  if (!type[1] %in% 1:5)
    stop("Unknown type sepecified!")
  
  s <- sum(moa@eig > tol)
  if (!is.null(n))
    s <- min(n, s)  
  ei <- (moa@eig/moa@eig[1])[1:s]
  m <- as.matrix(moa@ctr.tab[, 1:s])
  if (is.null(legend))
    legend <- rownames(m)
  

  if (type == 1) {

    barplot(moa@eig[1:s], col=col, ...)
    if (!is.null(axes))
      warning("type 1 renders barplot, axes argument are ignored")

  } else if( type == 2) { 

    barplot(m, legend.text=legend, plot=TRUE, col=col, ...)
    x <- barplot(m, plot=FALSE)
    lines(x, ei, pch=pch)
    if (!is.null(axes))
      warning("type 2 renders barplot, axes argument are ignored")

  } else if (type == 3) {

    barplot(m, beside=TRUE, legend.text=legend, col=col, ...)
    if (!is.null(axes))
      warning("type 3 renders barplot, axes argument are ignored")

  }  else if (type == 4) {

    if (is.null(col)) col <- 1
    matplot(t(m), type="l", col=col, lty=lty, ...)
    legend(x=lg.x, y=lg.y, legend=legend, col=col, lty=lty)
    if (!is.null(axes))
      warning("type 4 renders lines, axes argument are ignored")

  } else if (type == 5) {

    if (!inherits(axes, "numeric") | length(axes) != 2)
      stop("axes must consists of two integers")
    if (is.null(col)) col <- 1
    plot(m[, axes], col=col, pch=pch, ...)
    legend(x=lg.x, y=lg.y, legend=legend, pch=pch, col=col)
  }
}


# ==============================================================================
# ==                                                                          ==
# ==               RV coeff                                                   ==
# ==                                                                          ==
# ==============================================================================
.heatmap.rv <- function(moa, ...) {
  heatmap(moa@RV, symm=TRUE, ...)
}


