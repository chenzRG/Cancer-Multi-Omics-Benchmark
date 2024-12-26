\name{moa}
\alias{moa}

\title{
  Multiple omics data analysis using MFA or STATIS
}
\description{
    Analysis multiple omics data using MFA or STATIS. The input multiple tables are in a 
    form that columns are samples and rows are variables/features.
}
\usage{
  moa(data, proc.row="center_ssq1", w.data="inertia", w.row=NULL, statis=FALSE, moa=TRUE)
}

\arguments{
  \item{data}{
    A list of \code{data.frame} or \code{matrix} that contains the 
    input datas, the columns in all datasets should be samples/observations 
    (which need to be matched) and rows should be variables.
  }
  \item{proc.row}{
    Preprocessing of rows of datasets, should be one of 
    \code{none} - no preprocessing, 
    \code{center} - center only, 
    \code{center_ssq1} - center and scale (sum of squred values equals 1), 
    \code{center_ssqN} - center and scale (sum of squred values equals 
      the number of columns),
    \code{center_ssqNm1} - center and scale (sum of squred values equals 
      the number of columns - 1)
    MFA corresponds to "proc.row=center_ssq1" and 'w.data="lambda1"'
  }
  \item{w.data}{
    The weights of each separate dataset, should be one of 
    \code{uniform} - no weighting, 

    \code{lambda1} - weighted by the reverse of the first eigenvalue of 
      each individual dataset

    or \code{inertia} - weighted by the reverse of the total inertia. 
    See detail.
  }
  \item{w.row}{
    If it is not null, it should be a list of positive numerical vectors, 
    the length of which should be the same with the number of rows of 
    each dataset to indicated the weight of rows of datasets.
  }
  \item{statis}{
    A logical indicates whether STATIS method should be used. See details.
  }
  \item{moa}{
    Logical; whether the output should be converted to an object of class \code{\link{moa-class}}
  }
}

\details{
    Different methods employs different precessing of row and datasets. 
    For multipple factorial analysis (MFA), the rows of each dataset are 
    first centered and scaled, then each dataset is weighted by the 
    reverse of its first eigenvalue  (proc.row=center_ssq1, 
    w.data="lambda1"). This algorithm does not have a well defined criterion to be optimized (see reference). 

    If statis=TRUE, the statis algorithm will be used, that is, each 
    dataset will be further weighted so that datasets closer to the 
    overall structure will receive a higher weight.
}

\value{
    An object of class \code{\link{moa-class}}.
}

\references{
    Herve Abdi, Lynne J. Williams, Domininique Valentin and Mohammed Bennani-Dosse. STATIS and DISTATIS: optimum multitable principal component analysis and three way metric multidimensional scaling. WIREs Comput Stat 2012. Volume 4, Issue 2, pages 124-167
    Herve Abdi, Lynne J. Williams, Domininique Valentin. Multiple factor analysis: principal component analysis for multitable and multiblock data sets. WIREs Comput Stat 2013
}

\author{
    Chen Meng
}

\seealso{
  \code{\link{sup.moa}}, \code{\link{mogsa}}. More about plot see \code{\link{moa-class}}.
}
\examples{

  # library(mogsa)
  # loading data
  data(NCI60_4arrays)
  # run analysis
  ana <- moa(NCI60_4arrays, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
  # plot
  # plot eigen value
  plot(ana, value = "eig", type = 2)
  # plot the normalized (percentage) eigen value
  plot(ana, value = "tau", type = 2)
  # ploting the observations
  colcode <- as.factor(sapply(strsplit(colnames(NCI60_4arrays$agilent), split="\\\."), "[", 1))
  plot(ana, type = 1, value = "obs", col=colcode)
  plot(ana, type = 2, value = "obs", col=colcode, data.pch=1:4)
  # plot variables/features in each data sets
  plot(ana, value = "var", layout=matrix(1:4, 2, 2))
  # plot the RV coefficients for the data sets
  plot(ana, value = "RV")

  # to extract the components representing significant concordance structures between datasets
  bt <- bootMoa(moa = ana, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE, B = 20)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
\keyword{ PCA }
\keyword{ MVA }
\keyword{ MFA }
\keyword{ STATIS }