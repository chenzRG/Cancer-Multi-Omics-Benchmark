\name{displayClusters}
\alias{displayClusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot given similarity matrix by clusters
}
\description{
Visualize the clusters in given similarity matrix
}
\usage{
displayClusters(W, group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{
Similarity matrix
}
  \item{group}{
A vector containing the labels for each sample in W. 
}
}
\value{
Plots given similarity matrix with patients ordered to form clusters.
}
\author{
Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

## First, set all the parameters:
K = 20;			# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 10; 		# Number of Iterations, usually (10~20)

## Data1 is of size n x d_1, 
## where n is the number of patients, d_1 is the number of genes, 
## Data2 is of size n x d_2, 
## where n is the number of patients, d_2 is the number of methylation
data(Data1)
data(Data2)

## Here, the simulation data (SNFdata) has two data types. They are complementary to each other. 
## And two data types have the same number of points. 
## The first half data belongs to the first cluster; the rest belongs to the second cluster.
truelabel = c(matrix(1,100,1),matrix(2,100,1)); ## the ground truth of the simulated data

## Calculate distance matrices
## (here we calculate Euclidean Distance, you can use other distance, e.g,correlation)

## If the data are all continuous values, we recommend the users to perform 
## standard normalization before using SNF, 
## though it is optional depending on the data the users want to use.  
# Data1 = standardNormalization(Data1);
# Data2 = standardNormalization(Data2);


## Calculate the pair-wise distance; 
## If the data is continuous, we recommend to use the function "dist2" as follows 
Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));

## next, construct similarity graphs
W1 = affinityMatrix(Dist1, K, alpha)
W2 = affinityMatrix(Dist2, K, alpha)

## These similarity graphs have complementary information about clusters.
displayClusters(W1, truelabel);
displayClusters(W2, truelabel);

}
