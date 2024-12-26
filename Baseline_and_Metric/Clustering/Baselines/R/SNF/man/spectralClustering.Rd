\name{spectralClustering}
\alias{spectralClustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spectral Clustering
}
\description{
Perform the famous spectral clustering algorithms. There are three variants. The default one is the third type. 
}
\usage{
spectralClustering(affinity, K, type = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{affinity}{
Similarity matrix
}
  \item{K}{
Number of clusters
}
  \item{type}{
The variants of spectral clustering to use.
}
}
\value{
A vector consisting of cluster labels of each sample.
}
\author{
Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{

## First, set all the parameters:
K = 20;##number of neighbors, usually (10~30)
alpha = 0.5; ##hyperparameter, usually (0.3~0.8)
T = 20; ###Number of Iterations, usually (10~50)

## Data1 is of size n x d_1, 
## where n is the number of patients, d_1 is the number of genes, 
## Data2 is of size n x d_2, 
## where n is the number of patients, d_2 is the number of methylation
data(Data1)
data(Data2)

## Calculate distance matrices (here we calculate Euclidean Distance, 
## you can use other distance, e.g. correlation)
Dist1 = dist2(as.matrix(Data1),as.matrix(Data1))
Dist2 = dist2(as.matrix(Data2),as.matrix(Data2))

## Next, construct similarity graphs
W1 = affinityMatrix(Dist1, K, alpha)
W2 = affinityMatrix(Dist2, K, alpha)

# Next, we fuse all the graphs
# then the overall matrix can be computed by
W = SNF(list(W1,W2), K, T)

## With this unified graph W of size n x n, 
## you can do either spectral clustering or Kernel NMF. 
## If you need help with further clustering, please let us know. 

## You can display clusters in the data by the following function
## where C is the number of clusters.
C = 2

## You can get cluster labels for each data point by spectral clustering
labels = spectralClustering(W, C)
}
