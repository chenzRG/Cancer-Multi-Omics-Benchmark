\name{moa.sup-class}
\Rdversion{1.1}
\docType{class}
\alias{moa.sup-class}

\title{Class \code{"moa.sup"}}
\description{
	moa.sup class desc.
}
\section{Objects from the Class}{
	Objects can be created by calls of the form \code{new("moa.sup", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{sup}:}{Object of class \code{"list"}, the matrix of supplementary data. }
    \item{\code{coord.sep}:}{The projection of geneset infromation on each separate data. }
    \item{\code{coord.comb}:}{The projection of geneset infromation on total dataset.}
    \item{\code{score}:}{the gene set-sample pathway score}
    \item{\code{score.data}:}{the gene set-sample pathway score, data separate}
    \item{\code{score.pc}:}{the gene set-sample pathway score, PC separate}
    \item{\code{score.sep}:}{the gene set-sample pathway score, separate. }
    \item{\code{p.val}:}{the p value matrix have the same dimension with score matrix. }
    \item{\code{p.val.corrected}:}{the matrix of corrected p values. }
  }
}
\section{Methods}{
	There is no generic function for objects of "moa.sup", but have specific function, including:
	- decompose.gs.ind
	- box.gs.feature
	- plotGS
	- decompose.gs.group
}

\author{
 Chen Meng
}

\seealso{
	objects to See Also as \code{\link{decompose.gs.ind}}, \code{\link{box.gs.feature}}, \code{\link{plotGS}}, \code{\link{decompose.gs.group}}.
}
\examples{
	showClass("moa.sup")
	data(NCI60_4array_supdata)
	data(NCI60_4arrays)

	sapply(NCI60_4array_supdata, dim)
	ana <- moa(NCI60_4arrays, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
	plot(ana, value="eig")
	smoa <- sup.moa(ana, sup=NCI60_4array_supdata, nf=5)
}
\keyword{classes}
