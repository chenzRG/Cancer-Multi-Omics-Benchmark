\name{moa-class}
\Rdversion{1.1}
\docType{class}
\alias{moa-class}
\alias{plot,moa,missing-method}

\title{Class \code{"moa"}}
\description{
    moa class object
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("moa", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{eig}:}{eigen values} 
    \item{\code{tau}:}{The percentage of explained variance by each datasets sparately. } 
    \item{\code{partial.eig}:}{matrix, rows indicate the partial eigenvalues from each data.} 
    \item{\code{eig.vec}:}{a matrix, eigenvectors.} 
    \item{\code{loading}:}{the coordinate of variables/features.} 
    \item{\code{fac.scr}:}{factor score of observations.} 
    \item{\code{partial.fs}:}{partial factor score. } 
    \item{\code{ctr.obs}:}{ contribution of each observation to the total factor score. } 
    \item{\code{ctr.var}:}{ contribution of each variables to the total variance. } 
    \item{\code{ctr.tab}:}{contribution of each data to the total variance. } 
    \item{\code{RV}:}{ pairwise RV coefficients  } 
    \item{\code{w.row}:}{ weight of rows }
    \item{\code{w.data}:}{ weight of datasets } 
    \item{\code{data}:}{ the original input data } 
    \item{\code{tab.dim}:}{ the dimension of each input data } 
    \item{\code{call}:}{ call } 
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "moa", y = "missing")}: 
    Argument "value" sould be one of "eig", "tau", "obs", "var" and "RV"
    

    if value = "eig", the eigenvalue would be plotted as scree plot. The
    following arguments could be set:

      
      type=1 - The type of plot to show eigenvalues.
        (type=1: the eigenvalue are plotted;
        type=2: partial eigenvalue shown as concatenated bars;
        type=3: partial eigenvalue shown as bars side by side;
        type=4: matplot view of eigenvales, lty need to be set;
        type=5; the two dimensional plot of partial eigenvalues, axes and 
          pch need to be set in this case.) 

      axes=NULL - The axes selected to plot 

      n=NULL - Top n eigenvalues to be drawn 

      tol=1e-5 - The tolerance of eigenvalue, eigenvalues lower than 
        this value will not be shown. 

      legend=NULL - legend to put, a character string as calling legend 
        function 

      col=NULL - The color of partial eigenvalues from each data set 

      lty=1 - The line type used in the matplot, used when type =4 

      pch=NULL - the pch to draw 2D partial eigen plot, when type = 5 
        used 

      lg.x="topright" - The position of legend 

      lg.y=NULL - Poistion argument passed to function "legend" 

      ... - other arguemnts passed to functions 

    if value = "tau", the same with eig, but in the eigenvalues are scaled to 1 


    if value = "obs", the observation space will be shown, the following argument could be set:


      axes=1:2 - Which axes should be draw

      type=1 - Which type, see below (for type=1: the center points draw; type=2: the separate factor scores linked by lines; ... will be passed to function "points")

      data.pch=20 - the pch of dataset, if type=1, the first one is used

      col=1 - the color of observations, recycled used by data.frame

      label=FALSE - A logical indicates if labels should be shown

      lg.x="topright" - Position of legend 

      lg.y=NULL - Position of legend 

      xlim=NULL - The x limit 

      ylim=NULL - The y limit 

      label.cex=1 - the cex of text 

      ... 

    var - the separate gene view, layout can be specified 
    
    RV - the heatmap of RV coefficients 
    }

    \item{show}{\code{signature(x = "moa", y = "missing")}:
      show "moa" object
     }
	 }
}
\references{
    Herve Abdi, Lynne J. Williams, Domininique Valentin and Mohammed Bennani-Dosse. STATIS and DISTATIS: optimum multitable principal component analysis and three way metric multidimensional scaling. WIREs Comput Stat 2012. Volume 4, Issue 2, pages 124-167

    Herve Abdi, Lynne J. Williams, Domininique Valentin. Multiple factor analysis: principal component analysis for multitable and multiblock data sets. WIREs Comput Stat 2013
}

\author{
    Chen Meng
}

\examples{
    showClass("moa")
    # load("R/mogsa/data/NCI60_4arrays.rda")
    data(NCI60_4arrays)
    ana <- moa(NCI60_4arrays, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)

    plot(ana, value="eig")
    plot(ana, value="tau", type=2)
}
\keyword{classes}
