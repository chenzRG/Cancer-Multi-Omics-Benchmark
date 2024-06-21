\name{bootMoa}
\alias{bootMoa}

\title{
  Significant components in "moa" returned by function "moa".
}
\description{
    Using bootstrap method to extract the components representing significant concordance structures between datasets from "moa" (returned by function "moa").
}
\usage{
  bootMoa(moa, proc.row="center_ssq1", w.data="inertia", w.row=NULL, statis=FALSE,
    mc.cores=1, B = 100, replace=TRUE, resample=c("sample", "gene", "total"),
    plot=TRUE, log="y", tol = 1e-7)
}

\arguments{
  \item{moa}{
    An object of \code{\link{moa}} returned by \code{\link{moa}}.
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
  \item{mc.cores}{
    Integer; number of cores used in bootstrap. This value is passed to function \code{mclapply}
  }
  \item{B}{
    Integer; number of bootstrap
  }
  \item{replace}{
    Logical; sampling with or without replacement
  }
  \item{resample}{
    Could be one of "sample", "gene" or "total". "sample" and "gene" means sample-wise 
    and variable-wise resampling, repectively. "total" means total resampling.
  }
  \item{plot}{
    Logical; whether the result should be plotted.
  }
  \item{log}{
    Could be "x", "y" or "xy" for plot log axis.
  }
  \item{tol}{
    The minimum eigenvalues shown in the plot. 
  }
}

\details{
    set plot=TRUE to help selecting significant components.
}

\value{
    A matrix where columns are components and rows are variance of PCs from bootstrap samples.
}

\references{
    Herve Abdi, Lynne J. Williams, Domininique Valentin and Mohammed Bennani-Dosse. STATIS and DISTATIS: optimum multitable principal component analysis and three way metric multidimensional scaling. WIREs Comput Stat 2012. Volume 4, Issue 2, pages 124-167
    Herve Abdi, Lynne J. Williams, Domininique Valentin. Multiple factor analysis: principal component analysis for multitable and multiblock data sets. WIREs Comput Stat 2013
}

\author{
    Chen Meng
}

\seealso{
  \code{\link{moa}}, \code{\link{sup.moa}}, \code{\link{mogsa}}. More about plot see \code{\link{moa-class}}.
}
\examples{
  # see function moa
}