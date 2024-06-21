\name{NCI60_4arrays}
\alias{NCI60_4arrays}
\docType{data}
\title{
  Microarray gene expression profiles of the NCI 60 cell lines from 
  4 different platforms
}
\description{
  The 60 human tumour cell lines are derived 
  from patients with leukaemia, melanoma, lung, colon, central 
  nervous system, ovarian, renal, breast and prostate cancers. The cell line
  panel is widely used in anti-cancer drug screen. In this dataset,
  a subset of microarray gene expression of the NCI 60 cell lines from 
  four different platforms are combined in a list, which could be used as 
  input to \code{mcia} directly.
}
\usage{data(NCI60_4arrays)}
\format{
  The format is: List of 4 \code{data.frame}s
\itemize{
  \item{\$agilent:}{\code{data.frame} containing 300 rows and 60 columns. 
  300 gene expression log ratio measurements of the NCI60 cell lines, by Agilent
  platform.}
  
  \item{\$hgu133:}{\code{data.frame} containing 298 rows and 60 columns. 
  298 gene expression log ratio measurements of the NCI60 cell lines, by H-GU133
  platform.}
  
  \item{\$hgu133p2:}{\code{data.frame} containing 268 rows and 60 columns. 
  268 gene expression log ratio measurements of the NCI60 cell lines, by H-GU133
  plus 2.0 platform.}

  \item{\$hgu95:}{\code{data.frame} containing 288 rows and 60 columns. 
  288 gene expression log ratio measurements of the NCI60 cell lines,  by H-GU95
  platform.}
  }
}

\value{
    \code{NCI60_4arrays} will be loaded in your working space. 
}

\source{
 Cell Miner \url{http://discover.nci.nih.gov/cellminer/}
}
\references{
  Reinhold WC, Sunshine M, Liu H, Varma S, Kohn KW, Morris J, Doroshow J, Pommier Y
  CellMiner: A Web-Based Suite of Genomic and Pharmacologic Tools 
  to Explore Transcript and Drug Patterns in the NCI-60 Cell Line Set. 
  Cancer Research. 2012 Jul, 15;72(14):3499-511
}

\keyword{datasets}
\keyword{NCI-60}
\keyword{Microarray}
