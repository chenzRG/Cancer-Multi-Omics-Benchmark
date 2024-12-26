biSoftK <- function (x, maxiter, kp, kt, weight.p, weight.t, pos = FALSE, unit.pb = TRUE, unit.tb = FALSE) {
  
  if (length(kp) < length(x)) 
    kp <- rep(kp, length.out = length(x))
  if (length(kt) < length(x))
    kt <- rep(kt, length.out = length(x))
  
  regproj <- function(xb, t, kp, kt, wp, wt, unit.pb, unit.tb, pos) {
    pb <- t(xb) %*% t/c(t(t) %*% t)
    if (!unit.pb)
      pb <- normvec(pb)
    pb <- softK(pb, kp, w = wp, pos = pos)
    if (unit.pb)
      pb <- normvec(pb)
    tb <- xb %*% pb
    tb <- softK(tb, kt, w = wt, pos = pos)
    if (unit.tb)
      tb <- normvec(tb)
    list(tb = tb, pb = pb)
  }
  
  t <- fast.svd(do.call("cbind", x))$u[, 1, drop = FALSE]
  for (i in 1:maxiter) { 
    told <- t
    rp <- mapply(SIMPLIFY = FALSE, function(x, kp, kt, wp, wt, upb, utb, pos) {
      regproj(x, t, kp, kt, unit.pb = upb, unit.tb = unit.tb, wp = wp, wt = wt, pos = pos) ## stopped here
    }, x = x, kp = kp, kt = kt, wp = weight.p, wt = weight.t, pos = pos, upb = unit.pb, utb = unit.tb)
    
    tm <- sapply(rp, "[[", "tb")
    w <- t(tm) %*% t/c(t(t) %*% t)
    w <- w/sqrt(sum(w^2))
    t <- tm %*% w
    if (isTRUE(all.equal(c(t), c(told)))) 
      break
    if (i == maxiter) 
      cat("  Note: maximum number of iterations was reached, algrithm may not converge.\n")
  }
  res <- list(tb = lapply(rp, "[[", "tb"), pb = lapply(rp, "[[", "pb"), t = t, w = w)
  return(res)
}

# # examples
# data("NCI60_4arrays")
# library(svd)
# library(corpcor)
# source("R/mogsa/R/concordance.R")
# source("R/mogsa/R/normvec.R")
# source("R/mogsa/R/softK.R")
# source("R/mogsa/R/biSoftK.R")
# source("R/mogsa/R/mbpca2.R")

# d <- lapply(NCI60_4arrays[2:4], function(x) scale(t(x)))
# # 
# x1 <- biSoftK(d, maxiter = 1000, kp = Inf, kt = Inf,
#               weight.p = rep(1, length(d)),
#               weight.t = rep(1, length(d)),
#               pos = FALSE,
#               unit.pb = TRUE, unit.tb = TRUE)
# 
# x1s <- biSoftK(d, maxiter = 1000, kp = Inf, kt = Inf,
#               weight.p = 1,
#               weight.t = 1,
#               pos = FALSE,
#               unit.pb = TRUE, unit.tb = TRUE)
# 
# 
# identical(x1, x1s)
# barplot(c(x1$t), col= as.factor(substr(colnames(NCI60_4arrays$agilent), 1, 2)))



# 
# x2 <- biSoftK(d, maxiter = 1000, kp = 30, kt = Inf,
#               weight.p = rep(1, length(d)),
#               weight.t = rep(1, length(d)),
#               pos = FALSE,
#               unit.pb = TRUE, unit.tb = TRUE)
# barplot(c(x2$t), col= as.factor(substr(colnames(NCI60_4arrays$agilent), 1, 2)))
# 
# x3 <- biSoftK(d, maxiter = 1000, kp = 40, kt = Inf,
#               weight.p = rep(1, length(d)),
#               weight.t = rep(1, length(d)),
#               pos = TRUE,
#               unit.pb = TRUE, unit.tb = TRUE)
# barplot(c(x3$t), col= as.factor(substr(colnames(NCI60_4arrays$agilent), 1, 2)))
# 
# 
# x4 <- biSoftK(d, maxiter = 1000, kp = 30, kt = 6,
#               weight.p = rep(1, length(d)),
#               weight.t = rep(1, length(d)),
#               pos = FALSE,
#               unit.pb = TRUE, unit.tb = TRUE)
# barplot(c(x4$t), col= as.factor(substr(colnames(NCI60_4arrays$agilent), 1, 2)),
#         names.arg = colnames(NCI60_4arrays$agilent), las = 2)
# 
# 
# plot(x4$t, x4$tb$hgu133)
# plot(x4$t, x4$tb$hgu133p2)
# plot(x4$t, x4$tb$hgu95)
# heatmap(t(d$hgu133[, which(x4$pb$hgu133 != 0) ]))
# barplot(t(d$hgu133[, which.max(x4$pb$hgu133) ]), las = 2)
# barplot(t(d$hgu133[, which.max(x4$pb$hgu133) ]), las = 2)

