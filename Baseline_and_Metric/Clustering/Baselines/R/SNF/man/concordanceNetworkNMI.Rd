\name{concordanceNetworkNMI}
\alias{concordanceNetworkNMI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Concordance Network NMI calculation
}
\description{
Given a list of affinity matrices, Wall, the number of clusters, return a matrix containing the NMIs between cluster assignments made with spectral clustering.
}
\usage{
concordanceNetworkNMI(Wall, C)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Wall}{
List of matrices. Each element of the list is a square, symmetric matrix that shows affinities of the data points from a certain view.
}
  \item{C}{
Number of clusters
}
}
\value{
Returns an affinity matrix that represents the neighborhood graph of the data points.
}
\author{
Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

# How to use SNF with multiple views

# Load views into list "dataL"
data(dataL)
data(label)

# Set the other parameters
K = 20 # number of neighbours
alpha = 0.5 # hyperparameter in affinityMatrix
T = 20 # number of iterations of SNF
# Normalize the features in each of the views.
#dataL = lapply(dataL, standardNormalization)

# Calculate the distances for each view
distL = lapply(dataL, function(x) dist2(x, x))

# Construct the similarity graphs
affinityL = lapply(distL, function(x) affinityMatrix(x, K, alpha))

# an example of how to use concordanceNetworkNMI
Concordance_matrix = concordanceNetworkNMI(affinityL, 3);

## The output, Concordance_matrix,
## shows the concordance between the fused network and each individual network. 

}
