# moa - multiple omics data analysis 
# 
# input arguments
# data 
#   a list of data.frame or matrix that rows represent variables (genes)
#   and coloumns represent measurements/observations (samples, cell lines)
# method
#   either "statis" or "mfa"
# full.analysis
#   a logical indicate wether the full analysis should be performed
#   set as TRUE in default. In the bootstrapping preocedure, it set as
#   FALSE
# 
# detail
#   if full analysis is FALSE, returns
#     eig - eigen values
#     eig.vec - eigen vector (columns, in omics data, samples/observations)
#     loading - eigen vector for rows (in omics data, genes/molecules)
#     fac.scr - factor score
#     tab.dim - the dimension of each table
#     method - which method, "statis" or "mfa"
#     full.analysis - whether full analysis is performed
#   if full analysis is TRUE, returns
#     eig, eig.vec, loading, tab.dim, method, full.analysis, fac.scr -see above
#     tau - explained variance (percentage) by each eigenvector
#     fac.scr - factor.scores, calculated see abdi. 2012 statis ...
#     partial.fs - partial factor scores
#     ctr.obs - contribution of observations/samples
#     ctr.var - contribution of variables
#     ctr.tab - contribution of tables
#     partial.eig - partial eigen value

moa <- function(data, proc.row="center_ssq1", w.data="inertia", w.row=NULL, statis=FALSE, moa=TRUE) {
  
  kd <- data  
  data <- lapply(data, as.matrix)
  nRows <- sapply(data, nrow)

  # check 
  if (!is.null(w.row)) {
    wr <- do.call("c", w.row)
    if (any(wr <= 0))
      stop ("postive row weights are required.")
    if (! all.equal(sapply(w.row, length), nRows))
      stop("w.row should be a list of vector that contain the weight for each row in data! The length of 
           the element of w.row should equals the # of rows in each data.")
  }
  
  checkRowName <- sapply(data, function(x) is.null(rownames(x)))
  checkColName <- sapply(data, function(x) is.null(colnames(x)))
  checkTabName <- is.null(names(data))
  if (any(c(checkRowName, checkColName, checkTabName)))
    stop ("Table name or Colnames or Rownames are missing in all/one data!")
  nCols <- sapply(data, ncol)
  if (length(unique(nCols)) != 1) 
    stop ("Unequal number of columns! Columns need to be matched.")
  
  # preprocessing of row 
  nObs <- unique(nCols)
  nData <- length(data)
  data <- lapply(data, .proc.row, method = proc.row)
  
  # create weight
  wObs <- rep(1/nObs, nObs)
  wD <- .w.data(data, method = w.data, statis=statis)
  wData <- rep(sqrt(wD$w), nRows)
  if (!is.null(w.row))
    wData <- wData * sqrt(wr)

  # data concatenation, weighting and decomposition
  d1 <- .concateTabs(wD$x)
  Xt <- d1$tab * wData / sqrt(nObs)
  sing <- svd(Xt)
  
  if (!moa)
    return(sing)

  decom <- .read.svd(x=sing, data=wD$x, 
                      M=wObs, A=wData^2,
                      design=rep(names(data), nRows))

  decom$data <- kd
  decom$tab.dim <- data.frame(sapply(data, dim), row.names=c("row", "col"))
  decom$call <- match.call()
  
  decom$proc.row <- proc.row
  decom$w.data <- w.data
  decom$w.row <- w.row
  
  colnames(decom$partial.eig) <- paste("PC", 1:ncol(decom$partial.eig), sep="")
  result <- new("moa",
                eig = decom$eig,
                tau = decom$tau,
                partial.eig = decom$partial.eig,
                eig.vec = decom$eig.vec,
                fac.scr = decom$fac.scr,
                loading = decom$loading,
                partial.fs = decom$partial.fs,
                ctr.obs = decom$ctr.obs,
                ctr.var = decom$ctr.var,
                ctr.tab = decom$ctr.tab,
                RV = decom$RV,
                w.row = decom$alpha,
                w.data = wData^2,
                data = decom$data,
                tab.dim = decom$tab.dim,
                call = decom$call)
  return(result)
}

