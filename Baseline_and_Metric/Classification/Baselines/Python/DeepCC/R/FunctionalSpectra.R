#' Preprocess Gene List
#'
#' This function preprocess gene list for futhur process
#'
#' @param geneList a named vector containing the value of gene expression
#' @return a named vecter containing the value of gene expression
#' @examples
#' preprocessGeneList(geneList)
preprocessGeneList <- function(geneList) {
  geneList <- geneList[which((!is.na(geneList)) & (names(geneList)!="") & (!is.na(names(geneList))))]
  geneList <- tapply(t(geneList), names(geneList), max)
  geneList[order(geneList, decreasing = TRUE)]
}

#' Calculate Enrichment Score
#'
#' This function calculates enrichment score of a gene list on a specific gene set.
#'
#' @param geneList a named vector containing the values of gene expression
#' @param geneSet a vector containing genes to represent a gene set
#' @return a numeric indicating enrichment score
#' @useDynLib DeepCC
#' @import Rcpp
#' @export
#' @examples
#' calcEnrichmentScore(geneList, geneSet)
calcEnrichmentScore <- function(geneList, geneSet)
{
  calcEnrichmentScoreCPP((names(geneList) %in% geneSet), geneList, 1)
}

#' Generate Functional Spectra
#'
#' This function generates functional spectra for given gene expressin profiles.
#'
#' @param eps a data.frame containing gene expression profiles (each row presents one sample)
#' @param geneSet a List containing gene sets (default: MSigDB v6)
#' @param cores a integer indicating cpu cores used in parallel computing (default = all cores -2 )
#' @return a data.frame containing functional spectra
#' @seealso  \code{\link{getFunctionalSpectrum}} for a single expression profile.
#' @importFrom foreach foreach %dopar%
#' @export
#' @examples
#' getFunctionalSpectra(eps)
getFunctionalSpectra <- function(eps, geneSets = 'MSigDBv7', scale = T, cores = parallel::detectCores() - 2) {
  if (geneSets == 'MSigDBv5') {
    data(MSigDBv5)
    geneSets = MSigDBv5
  } else if (geneSets == 'MSigDBv6') {
    data(MSigDBv6)
    geneSets = MSigDBv6
  } else if(geneSets == 'MSigDBv7') {
    data(MSigDBv7)
    geneSets = MSigDBv7
  }

  if(scale) eps <- scale(eps, scale = FALSE)

  doParallel::registerDoParallel(cores)
  res <- foreach(idx = 1:nrow(eps), .combine = rbind) %dopar% {
    geneList <- preprocessGeneList(eps[idx, ])
    sapply(geneSets, function(x) calcEnrichmentScore(geneList, x))
  }
  rownames(res) <- rownames(eps)
  res
}

#' Generate Functional Spectrum
#'
#' This function generates functional spectrum for a single gene expression profile.
#'
#' @param expressionProfile a named numeric vector containing gene expression profile
#' @param geneSets a List containing gene sets (default: MSiDB v6)
#' @param refExp a character indicating cancer typer according to TCGA's indentifier, or a named vector reference expression
#' @param logChange a logical flag indicating whether the input data is already in log change form, e.g., for two color microarray, you should turn it on. (default: FALSE)
#' @param inverseRescale a logical flag indicating whether we rescale the reference to the scale of input data. If your single sample is microarray data and the reference is RNA-Seq, you should turn it on. (default: FALSE)
#' @param filter a numeric indicating the cutoff value of expression. (default: -3)
#' @return a numeric vector containing functional spectrum
#' @note You can generate the reference expression profile from your previous data or public data, which is the same(similiar) cancer type and platform.
#' In DeepCC we also prepared average expression profiles of each cancer types in TCGA project as references. To use them, just use the TCGA identifier (COADREAD, BRCA, OV, etc.) to indicate the cancer type.
#' If your single sample is microarray data, we strongly sugguest turn the parameter \code{inverseRescale} on, since TCGA is RNA-Seq, which has very small expression value for low expressed genes, compared with microarray.
#' @seealso \code{\link{getFunctionalSpectra}} for a batch of gene expression profiles.
#' @export
#' @examples
#' getFunctionalSpectrum(ep, refExp = "COADREAD")
getFunctionalSpectrum <- function(expressionProfile, geneSets = 'MSigDBv7', refExp = NULL, logChange = F, inverseRescale = F, filter = -3) {
  expressionProfile <- unlist(expressionProfile)
  if(!logChange) {
    if(is.null(refExp)) stop("Must have a reference expression profile!")
    if(is.character(refExp)) {
      if(!(refExp %in% rownames(ref_eps))) stop(paste(refExp, "is not a support identifier of cancer types!\n Please use one in :", paste(row.names(ref_eps), collapse = " ")))
      refExp <- ref_eps[refExp, ]
    }
    # filter low expressed genes
    refExp <- refExp[refExp > filter]

    common <- intersect(names(expressionProfile), names(refExp))
    if(!inverseRescale) {
      expressionProfile <- predict(lm(refExp[common] ~ expressionProfile[common])) -expressionProfile[common]
    } else {
      expressionProfile <- expressionProfile[common] - predict(lm(expressionProfile[common] ~ refExp[common]))
    }
  }
  geneList <- preprocessGeneList(expressionProfile)

  if (geneSets == 'MSigDBv5') {
    data(MSigDBv5)
    geneSets = MSiDBv5
  } else if (geneSets == 'MSigDBv6') {
    data(MSigDBv6)
    geneSets = MSigDBv6
  } else if (geneSets == 'MSigDBv7') {
    data(MSigDBv7)
    geneSets = MSigDBv7
  }

  res <- sapply(1:length(geneSets), function(idx) calcEnrichmentScore(geneList, geneSets[[idx]]))
  names(res) <- names(geneSets)
  res
}
