\name{prepSupMoa}
\alias{prepSupMoa}

\title{
    Prepare sumpplementary tables for projection by sup.moa or mogsa.
}

\description{
    Convert a list of gene set information to a set of sumpplementary tables that can be used
    as input of function "sup.moa" or "mogsa".
}

\usage{
    prepSupMoa(X, geneSets, minMatch = 10, maxMatch = 500)
}

\arguments{
  \item{X}{
    A matrix/data.frame or a list of matrix/data.frame or a list of character vector. If it is a list of matrix/data.frame, row names of matrix/data.frame will be used to create the projection matrix. Otherwise the charater vectors will used to create the supplementary matirx.
}
  \item{geneSets}{
    Gene sets list or an object of class "GeneSet" or "GeneSetCollection". 
    A gene set list could be returned by prepGraphite or prepMolsigDB.
}
  \item{minMatch}{
    The minimum match of geneset.
}
  \item{maxMatch}{
    The maximum match genesets.
}
}

\details{
    Details here
}

\value{
    A list of matrix could used as supplementary tables by "sup.moa" or "mogsa".
}

\author{
    Chen Meng
}

\seealso{
  See Also as \code{\link{prepGraphite}} and \code{\link{prepMsigDB}}.
}

\examples{
  library(graphite)
  data(NCI60_4arrays)
  kegg <- pathways(species = "hsapiens", "kegg")
  pw <- c("Purine metabolism", "MAPK signaling pathway")
  gss <- prepGraphite(db = kegg[pw], id="symbol")
  gss <- lapply(gss, function(x) sub("SYMBOL:", "", x))
  sup_data1 <- prepSupMoa(NCI60_4arrays, geneSets=gss)
  gene_list <- lapply(NCI60_4arrays, rownames)
  sup_data2 <- prepSupMoa(gene_list, geneSets=gss)
}

