## ----style, eval=TRUE, echo=FALSE, results='asis'--------------------------
BiocStyle::latex(bibstyle="unsrt")
# BiocStyle::latex()

## ----include=FALSE---------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----include=FALSE---------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE, cache=TRUE, message=FALSE, out.width=".55\\textwidth", echo=TRUE, fig.width=6, fig.height=6, fig.align="center", result="markup", hold=TRUE
)

## ----loadLib---------------------------------------------------------------
# loading package and gene expression data
library(mogsa)
data(NCI60_4arrays)

## ----dataDim---------------------------------------------------------------
sapply(NCI60_4arrays, dim) # check dimensions of expression data

## ----auxiVar---------------------------------------------------------------
tumorType <- sapply(strsplit(colnames(NCI60_4arrays$agilent), split="\\."), "[", 1)
colcode <- as.factor(tumorType)
levels(colcode) <- c("red", "green", "blue", "cyan", "orange", 
                     "gray25", "brown", "gray75", "pink")
colcode <- as.character(colcode)

## ----mbpca1, fig.cap="The variance associated with each latent variable. Colors distinguishes the contributions from different data sets."----
moa <- mbpca(NCI60_4arrays, ncomp = 10, k = "all", method = "globalScore", option = "lambda1", 
             center=TRUE, scale=FALSE, moa = TRUE, svd.solver = "fast", maxiter = 1000)
plot(moa, value="eig", type=2)

## ----boot, fig.cap="permutation test"--------------------------------------
r <- bootMbpca(moa, mc.cores = 1, B=20, replace = FALSE, resample = "sample")

## ----mpbca2----------------------------------------------------------------
moas <- mbpca(NCI60_4arrays, ncomp = 3, k = 0.1, method = "globalScore", option = "lambda1", 
              center=TRUE, scale=FALSE, moa = TRUE, svd.solver = "fast", maxiter = 1000)

## ----scoreCor--------------------------------------------------------------
scr <- moaScore(moa)
scrs <- moaScore(moas)
diag(cor(scr[, 1:3], scrs))

## ----plot1, fig.width=10, fig.height=6-------------------------------------
layout(matrix(1:2, 1, 2))
plot(scrs[, 1:2], col=colcode, pch=20)
legend("topright", legend = unique(tumorType), col=unique(colcode), pch=20)
plot(scrs[, 2:3], col=colcode, pch=20)

## ----gap, fig.cap="gap statistic plot"-------------------------------------
gap <- moGap(moas, K.max = 12, cluster = "hcl")
layout(matrix(1, 1, 1))
gap$nClust

## ----cluster, fig.width=6, fig.height=3------------------------------------
hcl <- hclust(dist(scrs))
cls <- cutree(hcl, k=4)
clsColor <- as.factor(cls)
levels(clsColor) <- c("red", "blue", "orange", "pink")
clsColor <- as.character((clsColor))

heatmap(t(scrs[hcl$order, ]), ColSideColors = colcode[hcl$order], Rowv = NA, Colv=NA)
heatmap(t(scrs[hcl$order, ]), ColSideColors = clsColor[hcl$order], Rowv = NA, Colv=NA)

## ----coef------------------------------------------------------------------
genes <- moaCoef(moas)
genes$nonZeroCoef$agilent.V1.neg

## ----sessionInfo, results = 'asis', eval = TRUE, echo = TRUE---------------
toLatex(sessionInfo())

