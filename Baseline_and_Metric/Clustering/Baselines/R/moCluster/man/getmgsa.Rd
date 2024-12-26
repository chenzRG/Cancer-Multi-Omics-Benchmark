\name{getmgsa}
\alias{getmgsa}
\title{
  get values in an object of class "mgsa".
}
\description{
	get values/slot in an object of class "mgsa". The "mgsa" consists of two S4 class objects,
	\code{\link{moa-class}} and \code{\link{moa.sup-class}}. This function could extract 
	values in these two components directly.
}
\usage{
	getmgsa(mgsa, value)
}

\arguments{
  \item{mgsa}{
  	An object of class \code{\link{mgsa-class}}.
  }
  \item{value}{
	The name of the value want to extract from "mgsa". See detail for options.
  }
}
\details{
	if value in c("call", "moa", "sup"), the function equal function \code{\link{slot}}.
	
	if value in c("eig", "tau", "partial.eig", "eig.vec", "loading", 
	"fac.scr", "partial.fs", "ctr.obs", "ctr.var", "ctr.tab", "RV"), the function extact
	corresponding value from \code{\link{moa-class}}. 
	
	if value in c("data", "coord.sep", "coord.comb", "score", 
    "score.data", "score.pc", "score.sep", "p.val"), the function extract value from 
     \code{\link{moa.sup-class}}.
}
\value{
	The function return the selected value in "mgsa".
}

\references{
	TBA
}
\author{
	Chen Meng
}


\examples{
	# library(mogsa)
	# loading gene expression data and supplementary data
	data(NCI60_4array_supdata)
	data(NCI60_4arrays)
	mgsa <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=9,
	              proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
	part.eig <- getmgsa(mgsa, "partial.eig")
	barplot(as.matrix(part.eig))
}
