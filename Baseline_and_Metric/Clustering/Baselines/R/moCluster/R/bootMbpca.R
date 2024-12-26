bootMbpca <-
function(moa, mc.cores=1, B = 100, replace=TRUE, resample=c("sample", "gene", "total"), log="y",
                      ncomp = NULL, method = NULL, maxiter=1000, svd.solver=c("svd", "fast.svd", "propack"), plot=TRUE) {
  data <- moa@data
  call <- moa@call
  
  if (is.null(ncomp)) 
    ncomp <- call$ncomp
  if (is.null(method)) {
    method <- call$method
    cat(paste("method is set to '", method, "'.\n", sep = ""))
  }
  if (!is.null(call$maxiter) && call$maxiter > 0 && is.integer(call$maxiter)) {
    maxiter <- call$maxiter
    cat(paste("maxiter is set to ", maxiter, ".\n", sep = ""))
  }  
  if (is.null(call$k))
    call$k <- "all"
  if (call$k != "all") {
    call$k <- "all"
    call$verbose <- FALSE
    moa <- eval(call)
  }
  
  ncomp <- min(c(ncomp), length(moa@eig))
  svd.solver <- match.arg(svd.solver)
  resample <- match.arg(resample)
  btk <- bootMbpcaK(data, B = B, mc.cores=mc.cores, replace = replace, resample=resample,
                    option = "uniform", center = FALSE, scale = FALSE,
                    ncomp = ncomp, k = "all", method = method, maxiter=maxiter, svd.solver=svd.solver)
  
  if (plot) {
    sc <- min(ncol(btk), length(moa@eig))
    isc <- 1:sc
    
    boxplot(rbind(moa@eig[isc], btk[, isc]), col=NA, border = NA, log=log)
    
    sds <- apply(btk, 2, sd)
    means <- apply(btk, 2, mean)
    points(isc,  means, pch=15, col="red")
    points(isc,  means+1.96*sds, col="red", pch="_")
    points(isc,  means-1.96*sds, col="red", pch="_")
    segments(isc, means+1.96*sds, isc, means-1.96*sds, col="red")
    
    lines(isc, moa@eig[isc], pch=20)
    points(isc, moa@eig[isc], pch=20)
  }
  
  colnames(btk) <- paste("PC", 1:ncol(btk), sep="")
  rownames(btk) <- paste("sample", 1:nrow(btk), sep="")
  btk
}
