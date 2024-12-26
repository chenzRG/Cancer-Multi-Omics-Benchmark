\name{mbpca}
\alias{mbpca}

\title{
  Extension of PCA to analyze multiple data sets
}
\description{
  Three approaches are supplied in this function, consensus PCA (CPCA), 
  generalized CCA (GCCA) and multiple co-inertia analsyis (MCIA).
}
\usage{
mbpca(x, ncomp, method, k = "all", center = TRUE, 
    scale = FALSE, option = "uniform", maxiter = 1000, 
    moa = TRUE, verbose = TRUE, svd.solver = c("svd", "fast.svd", "propack"),
    k.obs = "all", w = NA, w.obs = NA,
    unit.p = FALSE, unit.obs = FALSE, pos = FALSE)
}

\arguments{
  \item{x}{
  A \code{list} of \code{matrix} or \code{data.frame}, where rows are variables
  and columns are samples. The columns among the matrices need to be match but the 
  variables do not need to be.
}
  \item{ncomp}{
  An integer; the number of components to calculate. To calculate more components 
  requires longer computational time.
}
  \item{method}{
  A character string could be one of c("globalScore", "blockScore", "blockLoading").
  The "globalScore" approach equals consensus PCA; The "blockScore" approach equals 
  generalized canonical correlation analysis (GCCA); The "blockLoading" approach 
  equals multiple co-inertia anaysis (MCIA); 
}
  \item{k}{
  The absolute number (if k >= 1) or the proportion (if 0<k<1) of non-zero coefficients 
  for the variable loading vectors. It could be a single value or a vector has the same length as x so 
  the sparsity of individual matrix could be different.
}
  \item{center}{
  Logical; if the variables should be centered
}
  \item{scale}{
  Logical; if the variables should be scaled
}
  \item{option}{
  A charater string could be one of c("lambda1", "inertia", "uniform") to indicate
  how the different matrices should be normalized. If "lambda1", the matrix is 
  divided by its the first singular value, if "inertia", the matrix is divided by
  its total inertia (sum of square), if "uniform", none of them would be done.
}
  \item{maxiter}{
  Integer; Maximum number of iterations in the algorithm
}
  \item{moa}{
  Logical; whether the output should be converted to an object of class \code{\link{moa-class}}
}
  \item{verbose}{
  Logical; whether the process (# of PC) should be printed
}
  \item{svd.solver}{
  A charater string could be one of c("svd", "fast.svd", "propack"). The default "fast.svd " has
  a good compromise between the robustness and speed. "propack" is the fastest but may failed to 
  converge in practice.
}
  \item{k.obs}{
  The absolute number (if k >= 1) or the proportion (if 0<k<1) of non-zero coefficients 
  for the observations. Sparse factor scores for observation are used by sparse concordance analysis.  
  (New arguments from v1.12)  
}
  \item{w}{
  The weight of variables. It could be given in the following format: 1) NA or a numeric value: all variables have the same weight; 2) A vector of numeric values, the vecter has the same length as x: variables in each block shares the same weight; 3) A list of vector, each vector in the list has the same length as the number of row in the corresponding table/block, then each variable use a different weight. See detail how to select weight.  (New arguments from v1.12)  
}
  \item{w.obs}{
  The weight of observations, see w. (New arguments from v1.12)  
}
  \item{unit.p}{
  A logical value, whether the loading vectors (for variables) for each table/block should be unit length.  
}
  \item{unit.obs}{ A logical value, whether the score vectors (for observations) for each table/block should be unit length. (New arguments from v1.12)  
}
  \item{pos}{
  A logical value, whether only retain non-negative coefficients in loading and score vectors.  (New arguments from v1.12)  
}
}
\details{
  Select of weight for variables:
  In omics data, it is often true that low intensity variables suffers more noise. Therefore, 
  The variables with higher intensities are more reliable. If we consider this, we can use the total 
  sum intensity of a variable (or a tranform of it) as weight, the model would prefer to select high 
  intensity variables.
}
\value{
    An object of class \code{\link{moa-class}} (if \code{moa=TRUE}) or 
    an \code{list} object contains the following elements:

    \code{tb} - the block scores 

    \code{pb} - the block loadings

    \code{t} - the global scores

    \code{w} - the wegihts of block scores to construct the global scor
}
\references{
  For clustering problem: Meng et al. 2015 moCluster: Identifying Joint Patterns Across Multiple Omics Data
  Sets. Journal of proteome research.
}
\author{
  Chen Meng
}
\note{
  no note
}


\seealso{
  see \code{\link{moa}} for non-iterative algorithms for multi-block PCA.
}
\examples{

data("NCI60_4arrays")
tumorType <- sapply(strsplit(colnames(NCI60_4arrays$agilent), split="\\\."), "[", 1)
colcode <- as.factor(tumorType)
levels(colcode) <- c("red", "green", "blue", "cyan", "orange", 
                     "gray25", "brown", "gray75", "pink")
colcode <- as.character(colcode)



moa <- mbpca(NCI60_4arrays, ncomp = 10, k = "all", method = "globalScore", option = "lambda1", 
             center=TRUE, scale=FALSE)
plot(moa, value="eig", type=2)
r <- bootMbpca(moa, mc.cores = 1, B=6, replace = FALSE, resample = "sample")

moas <- mbpca(NCI60_4arrays, ncomp = 3, k = 0.1, method = "globalScore", option = "lambda1", 
              center=TRUE, scale=FALSE)


scr <- moaScore(moa)
scrs <- moaScore(moas)
diag(cor(scr[, 1:3], scrs))

layout(matrix(1:2, 1, 2))
plot(scrs[, 1:2], col=colcode, pch=20)
legend("topright", legend = unique(tumorType), col=unique(colcode), pch=20)
plot(scrs[, 2:3], col=colcode, pch=20)

gap <- moGap(moas, K.max = 12, cluster = "hcl")
gap$nClust


hcl <- hclust(dist(scrs))
cls <- cutree(hcl, k=4)
clsColor <- as.factor(cls)
levels(clsColor) <- c("red", "blue", "orange", "pink")
clsColor <- as.character((clsColor))

heatmap(t(scrs[hcl$order, ]), ColSideColors = colcode[hcl$order], Rowv = NA, Colv=NA)
heatmap(t(scrs[hcl$order, ]), ColSideColors = clsColor[hcl$order], Rowv = NA, Colv=NA)

genes <- moaCoef(moas)
genes$nonZeroCoef$agilent.V1.neg


}

\keyword{ multi-blcok PCA }
\keyword{ GCCA }
\keyword{ CPCA }
\keyword{ MCIA }
