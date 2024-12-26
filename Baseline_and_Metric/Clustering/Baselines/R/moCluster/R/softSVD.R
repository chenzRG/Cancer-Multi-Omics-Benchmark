#' penalized SVD
#' @param x the matrix to be decomposed
#' @param nf the number of components want to computed
#' @param kv the number of nonzero coefficients for the right regularized sigular vectors
#' @param ku the number of nonzero coefficients for the left regularized sigular vectors
#' @param wv the weight for columns of x
#' @param wu the weight for rows of of x
#' @param pos Logical value, if only retain non-negative values in the 
#' @param maxiter the maximum number of iteration
#' @param tol the tolerance of error
#' @param verbose print progress of calulation
#' @author Chen Meng
#' @details The algorithm use a generalized version of NIPALS algorithm to allow 
#' sparsity, force non-negative values on both left and rigth singular vectors. 
#' In addition, different weghted could be introduced, the columns/rows 
#' with bigger weights are more likely to be selected. 
#' 
#' @description SVD with sparsity, non-negative and weight constrants
#' @seealso \code{\link{svd}}
#' @return the same as svd, list of three components, d, u, v
#' @examples 
#' library(corpcor)
#' source("R/mogsa/R/normvec.R")
#' source("R/mogsa/R/softK.R")
#' 
#' #' a random matrix
#' x <- matrix(rnorm(50), 5, 10)
#' ss <- svd(x)
#' ss2 <- softSVD(x, nf = 5)
#' ss$d
#' ss2$d
#' all.equal(abs(ss$u), abs(ss2$u))
#' all.equal(abs(ss$v), abs(ss2$v))
#' 
#' #' NCI60 data
#' library(mogsa)
#' data("NCI60_4arrays")
#' #
#' d <- as.matrix(NCI60_4arrays$agilent)
#' 
#' ###' compare with svd
#' ss <- svd(d, nu = 5, nv = 5)
#' ss2 <- softSVD(d, nf = 5)
#' ss$d
#' ss2$d
#' all.equal(abs(ss$u), abs(ss2$u))
#' all.equal(abs(ss$v), abs(ss2$v))
#' 
#' #' normalize the data
#' dt <- scale(t(d), scale = FALSE)
#' pp <- softSVD(x = dt, nf = 5, kv = 30, ku = 6, maxiter = 100)
#' barplot(pp$u[, 1], col = as.factor(substr(colnames(d), 1, 2)))
#' barplot(pp$v[, 1])
#' 
#' #' change sparsity
#' dt <- scale(t(d), scale = FALSE)
#' pp <- softSVD(x = dt, nf = 5, kv = 30, ku = 9, maxiter = 1000, pos = TRUE)
#' i <- 2
#' barplot(pp$u[, i], col = as.factor(substr(colnames(d), 1, 2)), las=2, names.arg = rownames(dt))
#' barplot(pp$v[, i])
#' 
#' #' change sparsity
#' pp <- softSVD(x = dt, nf = 5, kv = Inf, ku = 9, maxiter = 1000, pos = TRUE)
#' i <- 1
#' barplot(pp$u[, i], col = as.factor(substr(colnames(d), 1, 2)), las=2, names.arg = rownames(dt))
#' barplot(pp$v[, i])
#' 
#' #' change sparsity
#' pp <- softSVD(x = dt, nf = 5, kv = 30, ku = Inf, maxiter = 1000, pos = FALSE)
#' i <- 1
#' barplot(pp$u[, i], col = as.factor(substr(colnames(d), 1, 2)), las=2, names.arg = rownames(dt))
#' barplot(pp$v[, i])
#' 
#' #' use weight
#' w <- rowSums(d - min(d))
#' w <- w + max(w) #' prefer to select high intensity genes
#' 
#' pw <- softSVD(x = dt, nf = 6, kv = 30, ku = Inf, wv = w, maxiter = 1000, pos = FALSE, verbose = TRUE)
#' i <- 6
#' barplot(pw$u[, i], col = as.factor(substr(colnames(d), 1, 2)), las=2, names.arg = rownames(dt))
#' barplot(pw$v[, i])
#' 
#' 
#' i1 <- apply(pp$v, 1, function(x) any(x!=0))
#' i2 <- apply(pw$v, 1, function(x) any(x!=0))
#' plot(pw$u[, 1], pp$u[, 1])
#' layout(matrix(1:2, 1, 2))
#' hist(d[i1 & !i2, ], main = "No weight")
#' hist(d[i2 & !i1, ], main = "Weight")

softSVD <- function(x, nf = 1, kv = Inf, ku = Inf, wv = 1, wu = 1, pos = FALSE, 
                    maxiter = 50, tol = sqrt(.Machine$double.eps), verbose = FALSE) {
  
  regproj <- function(x, ku, kv, wu, wv, pos, maxiter) {
    u <- fast.svd(x)$u[, 1, drop = FALSE]
    for (i in 1:maxiter) {
      v <- crossprod(x, u)
      v <- softK(v, kv, w = wv, pos = pos)
      v <- normvec(v)
      uold <- u
      u <- x %*% v
      u <- softK(u, ku, w = wu, pos = pos)
      u <- normvec(u)
      if (isTRUE(all.equal(u, uold, check.attrbutes = FALSE))) {
        r <- list(d = attr(u, "length"), u = u, v = v)
        break()
      }
      if (i == maxiter) {
        r <- list(d = attr(u, "length"), u = u, v = v)
        warning("maxiter reached")
      }
    }
    r
  }
  res <- list(d = rep(NA, nf), u = matrix(NA, nrow(x), nf), v = matrix(NA, ncol(x), nf))
  for (i in 1:nf) {
    if (verbose)
      cat(paste("calculating component", i, "...\n"))
    r <- regproj(x, kv = kv, ku = ku, wv = wv, wu = wu, pos = pos, maxiter = maxiter)
    res$d[i] <- r$d
    res$u[, i] <- r$u
    res$v[, i] <- r$v
    x <- x - tcrossprod((r$u * r$d), r$v)
  }
  ii <- which(res$d > tol)
  res$d <- res$d[ii]
  res$u <- res$u[, ii, drop=FALSE]
  res$v <- res$v[, ii, drop=FALSE]
  res
}



