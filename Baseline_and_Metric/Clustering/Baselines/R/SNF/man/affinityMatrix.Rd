\name{affinityMatrix}
\alias{affinityMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Affinity matrix calculation
}
\description{
Computes affinity matrix from a generic distance matrix
}
\usage{
affinityMatrix(Diff, K, sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Diff}{
Distance matrix
}
  \item{K}{
Number of nearest neighbors
}
  \item{sigma}{
Variance for local model
}
}
\value{
Returns an affinity matrix that represents the neighborhood graph of the data points.
}
\references{
B Wang, A Mezlini, F Demir, M Fiume, T Zu, M Brudno, B Haibe-Kains, A Goldenberg (2014) Similarity Network Fusion: a fast and effective method to aggregate multiple data types on a genome wide scale. Nature Methods. Online. Jan 26, 2014  
}
\author{
Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{

## First, set all the parameters:
K = 20; ##number of neighbors, must be greater than 1. usually (10~30)
alpha = 0.5; ##hyperparameter, usually (0.3~0.8)
T = 20; ###Number of Iterations, usually (10~50)

## Data1 is of size n x d_1, 
## where n is the number of patients, d_1 is the number of genes, 
## Data2 is of size n x d_2, 
## where n is the number of patients, d_2 is the number of methylation
data(Data1)
data(Data2)

## Calculate distance matrices(here we calculate Euclidean Distance, 
## you can use other distance, e.g. correlation)
Dist1 = dist2(as.matrix(Data1),as.matrix(Data1))
Dist2 = dist2(as.matrix(Data2),as.matrix(Data2))

## Next, construct similarity graphs
W1 = affinityMatrix(Dist1, K, alpha)
W2 = affinityMatrix(Dist2, K, alpha)

}
