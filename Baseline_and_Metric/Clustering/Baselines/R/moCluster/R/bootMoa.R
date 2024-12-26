bootMoa <- function(
  moa, proc.row="center_ssq1", w.data="inertia", w.row=NULL, statis=FALSE,
  mc.cores=1, B = 100, replace=TRUE, resample=c("sample", "gene", "total"),
  plot=TRUE, log="y", tol = 1e-7
) {
  if (plot) {
    if (max(moa@eig) < tol)
      stop("tolerance is greater than highest eigenvalue, decrease tol!")
  }
  data <- lapply(moa@data, as.matrix)
  call <- moa@call
  
  checkParam <- function(call, x, input, opt) {
    if (is.null(call[[x]]))
      return()
    
    if (is.character(call[[x]])){
      v <- try(
        match.arg(call[[x]], opt), 
        silent = TRUE
      )
      if (class(v) == "try-error" || v != input)
        warning(paste0("Double check parameter '", x, "', it may not consistant with moa call!"))
    } else if (is.logical(call[[x]]))
      if (call[[x]] != input)
        warning(paste0("Double check parameter '", x, "', it may not consistant with moa call!"))
  }
  
  checkParam(call, "proc.row", input = proc.row, opt = c("none", "center", "center_ssq1", "center_ssqN", "center_ssqNm1"))
  checkParam(call, "w.data", input = w.data, opt = c("uniform", "lambda1", "inertia"))
  checkParam(call, "statis", input = statis, opt = c(TRUE, FALSE))
  checkParam(call, "w.row", input = w.row)
  
  resample <- match.arg(resample, c("sample", "gene", "total"))
  resampleMoa <- function(d, replace, resample, proc.row, w.data, w.row, statis) {
    rsd <- switch(
      resample,
      "sample" = lapply(d, function(x) structure(
        x[, sample(1:ncol(x), replace = replace)], 
        dimnames=list(rownames(x), colnames(x)))),
      "gene" = lapply(d, function(x) structure(
        t(apply(x, 1, sample, replace=replace)), 
        dimnames=list(rownames(x), colnames(x)))),
      "total" = lapply(d, function(x) structure(
        apply(x, 2, sample, replace=replace), 
        dimnames=list(rownames(x), colnames(x))))
      )
    
    res <- moa(data = rsd, proc.row=proc.row, w.data=w.data, w.row=w.row, statis=statis, moa=FALSE)
    res$d^2
  }
  r <- mclapply(1:B, mc.cores = mc.cores, function(x) 
    resampleMoa(data, replace = replace, resample = resample, 
                proc.row=proc.row, w.data=w.data, w.row=w.row, statis=statis)
  )
  btk <- do.call("rbind", r)
  if (plot) {
    sc <- min(ncol(btk), length(moa@eig))
    isc <- intersect(which(moa@eig > tol), 1:sc)
    boxplot(rbind(moa@eig[isc], btk[, isc]), col=NA, border = NA, log=log)
    sds <- apply(btk[, isc], 2, sd)
    means <- apply(btk[, isc], 2, mean)
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
