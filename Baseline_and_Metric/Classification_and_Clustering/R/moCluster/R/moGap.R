moGap <-
function(x, K.max, B=100, cluster=c("kmeans", "hclust"), plot=TRUE,
                  dist.method = "euclidean", dist.diag = FALSE, dist.upper = FALSE, dist.p = 2,
                  hcl.method = "complete", hcl.members = NULL,
                  km.iter.max = 10, km.nstart = 10, 
                  km.algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"), km.trace=FALSE) {
  
  cluster <- cluster[1]
  sr <- moaScore(x)
  fhclust <- function(x, k,
                      dist.method = "euclidean", dist.diag = FALSE, dist.upper = FALSE, dist.p = 2,
                      hcl.method = "complete", hcl.members = NULL) {
    d <- dist(x, method=dist.method, diag=dist.diag, upper=dist.upper, p=dist.p)
    hl <- hclust(d, method=hcl.method, members=hcl.members)
    cls <- cutree(hl, k=k)
    list(cluster=cls)
  }
  
  if (pmatch(cluster, "kmeans", nomatch = 0))
    v <- clusGap(sr, FUNcluster = kmeans, K.max = K.max, B = B, 
                iter.max=km.iter.max, nstart=km.nstart, algorithm=km.algorithm, trace=km.trace) else
      if (pmatch(cluster, "hclust", nomatch = 0))
        v <- clusGap(sr, FUNcluster = fhclust, K.max = K.max, B = B, 
                     dist.method = dist.method, dist.diag = dist.diag, dist.upper = dist.upper, dist.p = dist.p,
                     hcl.method = hcl.method, hcl.members = hcl.members) else
                       stop("Unkown clustering algorithm specified!")
  if (plot)
    plot(v)
  n <- sapply(c("firstSEmax", "Tibs2001SEmax", "globalSEmax","firstmax", "globalmax"), 
              function(x) maxSE(v$Tab[, 3], v$Tab[, 4], method = x))
  v$nClust <- n
  v
}
