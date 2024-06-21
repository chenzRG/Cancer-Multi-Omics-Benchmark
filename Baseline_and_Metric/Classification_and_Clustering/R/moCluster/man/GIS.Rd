\name{GIS}
\alias{GIS}
\title{
 calculate gene influential scores of genes in a gene set.
}
\description{
	Calculate the gene influential score of individual feature to the overall variance of GS score. 
	Using a leave-one-out procedure (See detail).
}
\usage{
	GIS(x, geneSet, nf=NA, barcol=NA, topN=NA, plot=TRUE, Fvalue=FALSE, ff=NA, cor=FALSE)
}

\arguments{
  \item{x}{
  	An object of class \code{\link{mgsa-class}}.
  }
  \item{geneSet}{
	A charater string or number to indicated the gene sets under 
	conserderation.
  }
  \item{nf}{
  	The number of PCs used in the caluculation of gene set scores. 
  	The default is NA, which means using all the PCs in the mogsa. This
  	should work for most of the cases.
  }
  \item{barcol}{
  	The color of the bars, which is used to distinguish features/genes from
  	different datasets, so its length should be the same as the number of
  	data sets.
  }
  \item{topN}{
	An positive integer specify the number of top influencers that should 
	to returned.
  }
  \item{plot}{
  	A logical indicate if the result should be plotted.
  }
  \item{Fvalue}{
  	A logical indicate if the GIS should be calculated in a supervised manner.
  }
  \item{ff}{
	The vector indicates the group of columns for calculating the F-ratio when
	Fvalue=TRUE. 
  }
  \item{cor}{
	A logical indicates whether use correlation between reconstructed expression with GSS. 
	This is faster than the standard GIS.
  }
}
\details{
	The evaluation of the importance of a single feature is calculated in 
	the supervised or unsupervised manner. 
	
	In the unsupervise manner, the value is calculated by:

	log2(var(GS_{-i})/var(GS))
	
	where GS is the gene set score, and the GS_{-i} is a recalculate of 
	gene set score without i'th feature. var() is the variance.

	In the supervised manner, the value is caluclated as the F-ratio over
	a class vector:

	log2(F(GS_{-i})/F(GS))

	Where F() is the calculation of F-ratio. The unsupervised GIS is encouraged
	since it works better for most of the cases in practice.
}
\value{
	An object of class \code{data.frame} contains three columns. The first column is the feature name, 
	the second columns is the gene influential score. The third columns indicates from where the
	feature/gene is selected. 
}

\references{
	TBA
}
\author{
	Chen Meng
}

\seealso{
	see \code{\link{annotate.gs}}
}
\examples{
	# library(mogsa)
	# loading gene expression data and supplementary data
	data(NCI60_4array_supdata)
	data(NCI60_4arrays)
	mgsa <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=9,
	              proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
	allgs <- colnames(NCI60_4array_supdata[[1]])

	# unsupervised measurement
	GIS(mgsa, allgs[1], topN = 5)

	# supervised measurement
	tissueType <- as.factor(sapply(strsplit(colnames(NCI60_4arrays$agilent), split="\\\."), "[", 1))
	GIS(mgsa, allgs[1], topN = 5, Fvalue = TRUE, ff = tissueType)
	# more PCs to calcualte
	GIS(mgsa, allgs[1], nf = 20, topN = 5, Fvalue = TRUE, ff = tissueType)
}
