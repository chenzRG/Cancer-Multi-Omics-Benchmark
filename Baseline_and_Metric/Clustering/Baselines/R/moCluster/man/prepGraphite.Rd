\name{prepGraphite}
\alias{prepGraphite}

\title{
    Prepare pathway gene sets from graphite package
}
\description{
    Prepare pathway gene sets from "graphite" package, which could be passed to "prepSupMoa" function. 
}
\usage{
    prepGraphite(db, id = c("entrez", "symbol"))
}
\arguments{
  \item{db}{
    The database to be used, an object of class either 'PathwayList' create by "pathways" function.
}
  \item{id}{
    Which identifier for output, either "entrez" or "symbol".
}
}
\details{
   Only support "entrez" or "symbol" output currently.
}
\value{
    This function returns an object of list containing gene set information, which could be further processed by function "prepSupMoa" to convert to the object that can be used as input of "sup.moa" or "mogsa".
}
\references{
    Sales G, Calura E and Romualdi C (2014). graphite: GRAPH Interaction from pathway Topological Environment. R package version 1.10.1. 
}
\author{
    Chen Meng
}

\seealso{
 See Also as \code{\link{prepMsigDB}} and \code{\link{prepSupMoa}}.
}
\examples{
  library(graphite)
  keggdb <- prepGraphite(db = pathways("hsapiens", "kegg")[1:3], id = "entrez")
}

\keyword{ pahtways }
\keyword{ graphite }
