\name{prepMsigDB}
\alias{prepMsigDB}

\title{
	Conver gmt format file to a list
}
\description{
	Convert a gmt file (Could be downloaded from MSigDB) to a list of gene sets information.
}
\usage{
	prepMsigDB(file)
}

\arguments{
  \item{file}{
  	The directory and file name of the gmt file.
}
}

\value{
	This function returns an object of list containing gene set information, which could be further processed by function "prepSupMoa" to convert to the object that can be used as input of "sup.moa" or "mogsa".
}

\author{
	Chen Meng
}

\seealso{
 See Also as \code{\link{prepGraphite}} and \code{\link{prepSupMoa}}.
}

\examples{
	# not run
	dir <- system.file(package = "mogsa")
	preGS <- prepMsigDB(file=paste(dir, 
		"/extdata/example_msigdb_data.gmt.gz", sep = ""))
}

