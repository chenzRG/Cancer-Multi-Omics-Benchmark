softK <- function(x, k, w = 1, pos = FALSE) {
  if (k <= 0)
    stop("k should be postive integers or (0, 1)")
  n <- length(x)
  if (k >= n) {
    r <- x
  } else {
    if (k < 1)
      k <- round(n * k)
    if (pos) {
      if (sum(x > 0) <= k) {
        r <- x
        r[r < 0] <- 0
        return(r)
      }
      ax <- x
    } else
      ax <- abs(x)
    
    sx <- sign(x)
    naxw <- -ax*w
    rk <- rank(naxw, ties.method = "min")
    i <- rk <= k
    maxi <- min(ax[i])
    off <- max(0, max(ax[ax < maxi]))
    r <- ax - off
    r[!i] <- 0
    r <- sx*r
  }
  return(r)
}

# # unit test
# library(RUnit)
# softK.ui <- function() {
#   checkEquals(softK(x = 1:7, k = 3), c(0, 0, 0, 0, 1, 2, 3))
#   checkEquals(softK(x = 1:7, k = 3), c(0, 0, 0, 0, 1, 2, 3))
#   checkEquals(softK(x = 1:6, k = 2, w = c(3, 3, 3, 1, 1, 1)), c(0, 1, 2, 0, 0, 5))
#   checkEquals(softK(x = 1:5, k = 2, w = c(3, 3, 3, 1, 1)), c(0, 1, 2, 0, 0))
#   checkEquals(softK(x = -3:5, k = 5), c(-2, -1, 0, 0, 0, 1, 2, 3, 4))
#   checkEquals(softK(x = -3:5, k = 6), c(-2, -1, 0, 0, 0, 1, 2, 3, 4))
#   checkEquals(softK(x = -3:5, k = 7), -3:5)
#   checkEquals(softK(x = -3:5, k = 8), -3:5)
#   checkEquals(softK(x = -3:3, k = 3, w = c(1, 3, 1, 1, 5, 1, 1)), c(-3, -2,  0,  0,  1, 0, 3))
#   checkEquals(softK(x = -3:3, k = 3, w = c(1, 3, 1, 1, 5, 1, 1), pos = TRUE), c(0, 0, 0, 0, 1, 2, 3))
#   checkEquals(softK(x = -3:3, k = 2, w = c(1, 3, 1, 1, 5, 1, 1), pos = TRUE), c(0, 0, 0, 0, 1, 0, 3))
#   checkEquals(softK(x = -3:3, k = 1, w = c(1, 3, 1, 1, 5, 1, 1), pos = TRUE), c(0, 0, 0, 0, 1, 0, 0))
#   checkEquals(softK(x = -3:3, k = 1, pos = TRUE), c(0, 0, 0, 0, 0, 0, 1))
#   checkEquals(softK(x = -3:3, k = 2, pos = TRUE), c(0, 0, 0, 0, 0, 1, 2))
#   checkEquals(softK(x = -3:3, k = 3, pos = TRUE), c(0, 0, 0, 0, 1, 2, 3))
#   checkEquals(softK(x = -3:3, k = 4, pos = TRUE), c(0, 0, 0, 0, 1, 2, 3))
# 
# }
# softK.ui()

