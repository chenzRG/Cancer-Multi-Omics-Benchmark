

########################################
###############   NEMO   ###############
########################################

raw.data.to.similarity <- function(raw.data) {
  return(lapply(raw.data, function(x) {return(cor(as.matrix(x), as.matrix(x)))}))
}

normalize.matrix <- function(data.matrix) {
  temp = data.matrix - rowMeans(data.matrix)
  should.keep = (apply(temp, 1, sd) != 0)
  return ((temp / apply(temp, 1, sd))[should.keep, ])
}

#' @title NEMO num clusters
#' @name nemo.num.clusters
#' @description Estimates the number of clusters in an affinity graph.
#' @param W the affinity graph.
#' @param NUMC possible values for the number of clusters. Defaults to 2:15.
#' @return the estimated number of clusters in the graph.
#' @export
nemo.num.clusters <- function(W, NUMC=2:15) {
  if (min(NUMC) == 1) {
    warning("Note that we always assume there are more than one cluster.")
    NUMC = NUMC[NUMC > 1]
  }
  W = (W + t(W))/2
  diag(W) = 0
  if (length(NUMC) > 0) {
    degs = rowSums(W)
    degs[degs == 0] = .Machine$double.eps
    D = diag(degs)
    L = D - W
    Di = diag(1/sqrt(degs))
    L = Di %*% L %*% Di
    print(dim(L))
    eigs = eigen(L)
    eigs_order = sort(eigs$values, index.return = T)$ix
    eigs$values = eigs$values[eigs_order]
    eigs$vectors = eigs$vectors[, eigs_order]
    eigengap = abs(diff(eigs$values))
    eigengap = (1:length(eigengap)) * eigengap

    t1 <- sort(eigengap[NUMC], decreasing = TRUE, index.return = T)$ix
    return(NUMC[t1[1]])
  }
}

#' @title Spectral clustering
#' @name spectralClustering
#' @export
spectralClustering = SNFtool::spectralClustering


#' @title Num neighbors ratio
#' @name NUM.NEIGHBORS.RATIO
#' @export
NUM.NEIGHBORS.RATIO = 6

#' @title NEMO affinity graph
#' @name nemo.affinity.graph
#' @description Constructs a single affinity graph measuring similarity across different omics.
#' @param raw.data A list of the data to be clustered, where each an entry is a matrix of features x samples.
#' @param k The number of neighbors to use for each omic. It can either be a number, a list of numbers
#' or NA. If it is a number, this is the number of neighbors used for all omics. If this is a list,
#' the number of neighbors are taken for each omic from that list. If it is NA, each omic chooses the
#' number of neighbors to be the number of samples divided by NUM.NEIGHBORS.RATIO.
#' @return A single matrix measuring similarity between the samples across all omics.
#' @export
nemo.affinity.graph <- function(raw.data, k=NA) {
  if (is.na(k)) {
    k = as.numeric(lapply(1:length(raw.data), function(i) round(ncol(raw.data[[i]]) / NUM.NEIGHBORS.RATIO)))
  } else if (length(k) == 1) {
    k = rep(k, length(raw.data))
  }
  sim.data = lapply(1:length(raw.data), function(i) {affinityMatrix(dist2(as.matrix(t(raw.data[[i]])),
                                                                as.matrix(t(raw.data[[i]]))), k[i], 0.5)})
  affinity.per.omic = lapply(1:length(raw.data), function(i) {
    sim.datum = sim.data[[i]]
    non.sym.knn = apply(sim.datum, 1, function(sim.row) {
      returned.row = sim.row
      threshold = sort(sim.row, decreasing = T)[k[i]]
      returned.row[sim.row < threshold] = 0
      row.sum = sum(returned.row)
      returned.row[sim.row >= threshold] = returned.row[sim.row >= threshold] / row.sum
      return(returned.row)
    })
    sym.knn = non.sym.knn + t(non.sym.knn)
    return(sym.knn)
  })
  patient.names = Reduce(union, lapply(raw.data, colnames))
  num.patients = length(patient.names)
  returned.affinity.matrix = matrix(0, ncol = num.patients, nrow=num.patients)
  rownames(returned.affinity.matrix) = patient.names
  colnames(returned.affinity.matrix) = patient.names

  shared.omic.count = matrix(0, ncol = num.patients, nrow=num.patients)
  rownames(shared.omic.count) = patient.names
  colnames(shared.omic.count) = patient.names

  for (j in 1:length(raw.data)) {
    curr.omic.patients = colnames(raw.data[[j]])
    returned.affinity.matrix[curr.omic.patients, curr.omic.patients] = returned.affinity.matrix[curr.omic.patients, curr.omic.patients] + affinity.per.omic[[j]][curr.omic.patients, curr.omic.patients]
    shared.omic.count[curr.omic.patients, curr.omic.patients] = shared.omic.count[curr.omic.patients, curr.omic.patients] + 1
  }

  final.ret = returned.affinity.matrix / shared.omic.count
  lower.tri.ret = final.ret[lower.tri(final.ret)]
  final.ret[shared.omic.count == 0] = mean(lower.tri.ret[!is.na(lower.tri.ret)])

  return(final.ret)
}

#' @title NEMO clustering
#' @name nemo.clustering
#' @description Performs multi-omic clustering on a datset using the NEMO algorithm.
#' Uses nemo.num.clusters to estimate the number of clusters.
#' @param omics.list A list of the data to be clustered, where each an entry is a matrix of features x samples.
#' @param k The number of neighbors to use for each omic. It can either be a number, a list of numbers
#' or NA. If it is a number, this is the number of neighbors used for all omics. If this is a list,
#' the number of neighbors are taken for each omic from that list. If it is NA, each omic chooses the
#' number of neighbors to be the number of samples divided by NUM.NEIGHBORS.RATIO.
#' @return A single matrix measuring similarity between the samples across all omics.
#' @export
nemo.clustering <- function(omics.list, num.clusters=NULL, num.neighbors=NA) {
  if (is.null(num.clusters)) {
    num.clusters = NA
  }

  graph = nemo.affinity.graph(omics.list, k = num.neighbors)
  if (is.na(num.clusters)) {
    num.clusters = nemo.num.clusters(graph)
  }
  clustering = spectralClustering(graph, num.clusters)
  names(clustering) = colnames(graph)
  return(clustering)
}
