bootMbpcaK <-
function(data, replace, B=100, mc.cores=1, resample = c("sample", "total", "gene"), 
                       ncomp, method, k, 
                       center=FALSE, scale=FALSE, option="uniform", 
                       maxiter=1000, svd.solver=c("svd", "fast.svd", "propack")) {
  
  
  resampleMbpca <- function(d, ncomp, method, k, center, scale, option, 
                            maxiter, replace, resample, svd.solver) {
    rsd <- switch(resample,
                  "sample" = lapply(d, function(x) x[, sample(1:ncol(x), replace = replace)]),
                  "gene" = lapply(d, function(x) t(apply(x, 1, sample, replace=replace))),
                  "total" = lapply(d, function(x) apply(x, 2, sample, replace=replace)))
                  # "total" = lapply(d, function(x) x[sample(1:nrow(x), replace = replace), sample(1:ncol(x), replace = replace)]))
    res <- mbpca(x = rsd, verbose = FALSE, moa=FALSE, 
                 ncomp=ncomp, method=method, k=k, center=center, 
                 scale=scale, option=option, maxiter=maxiter, svd.solver)
    diag(crossprod(res$t))
  }
  
  svd.solver <- match.arg(svd.solver)
  resample <- match.arg(resample)
  r <- mclapply(1:B, mc.cores = mc.cores, function(x) 
    resampleMbpca(data, ncomp, method, k, center, scale, option, maxiter, replace, resample, svd.solver))
  do.call("rbind", r)
  
}
