\name{NCI60_4array_supdata}
\alias{NCI60_4array_supdata}
\docType{data}
\title{
  supp data for Microarray gene expression profiles of the NCI 60 cell lines from 
  4 different platforms
}
\description{
  Supplmentary to NCI60_4arrays.
}
\usage{data(NCI60_4arrays)}
\format{
  The format is: List of 4 \code{matrix}
\itemize{
  \item{\$agilent:}{\code{matrix} containing 300 rows and 60 columns. 
  300 gene expression log ratio measurements of the NCI60 cell lines, by Agilent
  platform.}
  
  \item{\$hgu133:}{\code{matrix} containing 298 rows and 60 columns. 
  298 gene expression log ratio measurements of the NCI60 cell lines, by H-GU133
  platform.}
  
  \item{\$hgu133p2:}{\code{matrix} containing 268 rows and 60 columns. 
  268 gene expression log ratio measurements of the NCI60 cell lines, by H-GU133
  plus 2.0 platform.}

  \item{\$hgu95:}{\code{matrix} containing 288 rows and 60 columns. 
  288 gene expression log ratio measurements of the NCI60 cell lines,  by H-GU95
  platform.}
  }
}

\value{
    \code{NCI60_4array_supdata} will be loaded in your working space. 
}

\keyword{datasets}
\keyword{NCI-60}
\keyword{Microarray}
