.proc.row <- function (x, method="none") {
  if (! method %in% c("none", "center" ,"center_ssq1", "center_ssqN", "center_ssqNm1"))
    stop ("row should be in c(\"none\", \"center\", \"center_ssq1\", \"center_ssqN\", \"center_ssqNm1\")")
  
  x <- as.matrix(x)
  n <- ncol(x)
  
  center_ssq1 <- function(x) {
    x <- x-rowMeans(x)
    x <- x/sqrt(rowSums(x^2))
    x[is.nan(x)] <- 0
    return(x)
  }
  
  center_ssqn <- function(x) {
    x <- x-rowMeans(x)
    x <- x/sqrt(rowSums(x^2))
    x <- x * sqrt(n)
    x[is.nan(x)] <- 0
    return(x)
  }
  
  center_ssqnm1 <- function(x) {
    x <- x-rowMeans(x)
    x <- x/sqrt(rowSums(x^2))
    x <- x * sqrt(n-1)
    x[is.nan(x)] <- 0
    return(x)
  }
  
  x <- switch(method, 
              none = x,
              center = x-rowMeans(x),
              center_ssq1 = center_ssq1(x),
              center_ssqN = center_ssqn(x), 
              center_ssqNm1 = center_ssqnm1(x))
  return (x)
  
}

.w.data <- function (x, method="uniform", statis=FALSE) {
  # x is a list of data.frame/matrix
  
  if (! method %in% c("uniform", "lambda1", "inertia"))
    stop("unknow weighting method for data!")
  
  x.ori <- x
  nData <- length(x)
  
  w <- switch(method[1],
              "uniform" = rep(1, nData),
              "lambda1" = 1/sapply(x, function (d) eigen(crossprod(d))$values[1]),
              "inertia" = 1/sapply(x, function (d) sqrt(sum(d^2))))
  alpha <- w

  if (statis) {
    x <- mapply(SIMPLIFY = FALSE, function (x, w) {
      x * w
    }, x=x, w=w)

    ss <- lapply(x, crossprod)
    cc <- sapply(ss, c)
    cc <- crossprod(cc)
    svdc <- svd(cc)
    alpha <- svdc$u[, 1] / sum(svdc$u[, 1])
    alpha <- w^2 * alpha
  }
  list (w=alpha, x=x.ori)
}

.concateTabs <- function(tabs) {
  design <- sapply(tabs, nrow)
  names(design) <- names(tabs)
  tab <- do.call("rbind", tabs)
  list(tab=tab, design=design)
}


.read.svd <- function(x, data, M, A, design) {
  res <- list()
  res$eig <- x$d^2
  res$eig.vec <- x$v / sqrt(M)
  rownames(res$eig.vec) <-  colnames(data[[1]])
  
  res$loading <- as.data.frame(x$u / sqrt(A))
  rownames(res$loading) <- paste(unlist(lapply(data, rownames)), design, sep="_")
  res$fac.scr <- sweep(res$eig.vec, 2, x$d, "*")
  rownames(res$fac.scr) <- colnames(data[[1]])
  colnames(res$fac.scr) <- paste("PC", 1:ncol(res$fac.scr), sep = "")
  res$RV <- pairwise.rv(data, match="col")
  res$tau <- res$eig/sum(res$eig)
  Ql <- split(as.data.frame(res$loading), design)
  Ql <- Ql[names(data)]
  alpha <- split(A, design)
  alpha <- alpha[names(data)]
  
  res$partial.fs <- mapply(function(x, y, a) {
    t(a * x) %*% as.matrix(y)
  }, x=data, y=Ql, a=alpha, SIMPLIFY=FALSE)
  
  res$partial.fs <- lapply(res$partial.fs, as.data.frame)
  
  res$ctr.obs <- data.frame(sweep(M * res$fac.scr^2, 2, res$eig, "/"))
  res$ctr.var <- data.frame(A * res$loading^2)
  res$ctr.tab <- by(res$ctr.var, design, colSums)
  res$alpha <- alpha
  res$ctr.tab <- data.frame(t(sapply(res$ctr.tab, function(x) 
    x)), row.names=names(data))
  res$partial.eig <- data.frame(sweep(res$ctr.tab, 2, res$eig, "*"), 
                                row.names=names(data))
  return(res)
}