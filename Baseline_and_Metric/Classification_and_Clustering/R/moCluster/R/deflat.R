deflat <-
function(x, t, tb, pb, method="globalScore") {
  # globalScore, blockScore, blockLoading
  switch(method,
         "globalScore" = lapply(x, function(xb) { xb - t %*% t(t) %*% xb / c(t(t) %*% t) }),
         "blockLoading" = mapply(SIMPLIFY = FALSE, function(xb, pb) {xb - xb %*% pb %*% t(pb)}, x, pb),
         "blockScore" = mapply(SIMPLIFY = FALSE, function(xb, tb) {xb - tb %*% t(tb) %*% xb / c(t(tb) %*% tb)}, x, tb))
}
