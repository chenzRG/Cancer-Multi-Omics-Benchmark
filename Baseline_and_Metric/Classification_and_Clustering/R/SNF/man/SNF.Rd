\name{SNF}
\alias{SNF}
\title{
Similarity Network Fusion
}
\description{
Similarity Network Fusion takes multiple views of a network and fuses them together to construct an overall status matrix. The input to our algorithm can be feature vectors, pairwise distances, or pairwise similarities. The learned status matrix can then be used for retrieval, clustering, and classification.
}
\usage{
SNF(Wall, K, t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Wall}{
List of matrices. Each element of the list is a square, symmetric matrix that shows affinities of the data points from a certain view.
}
  \item{K}{
Number of neighbors in K-nearest neighbors part of the algorithm.
}
  \item{t}{
Number of iterations for the diffusion process.
}
}
\value{
W is the overall status matrix derived
}
\references{
B Wang, A Mezlini, F Demir, M Fiume, T Zu, M Brudno, B Haibe-Kains, A Goldenberg (2014) Similarity Network Fusion: a fast and effective method to aggregate multiple data types on a genome wide scale. Nature Methods. Online. Jan 26, 2014  

Concise description can be found here: 
http://compbio.cs.toronto.edu/SNF/SNF/Software.html
}
\author{
Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

## First, set all the parameters:
K = 20;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 20; 	# Number of Iterations, usually (10~20)

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
displayClusters(W1,truelabel);
displayClusters(W2,truelabel);

## next, we fuse all the graphs
## then the overall matrix can be computed by similarity network fusion(SNF):
W = SNF(list(W1,W2), K, T)

## With this unified graph W of size n x n, 
## you can do either spectral clustering or Kernel NMF. 
## If you need help with further clustering, please let us know. 

## You can display clusters in the data by the following function
## where C is the number of clusters.
C = 2 								# number of clusters
group = spectralClustering(W,C); 	# the final subtypes information
displayClusters(W, group)

## You can get cluster labels for each data point by spectral clustering
labels = spectralClustering(W, C)

plot(Data1, col=labels, main='Data type 1')
plot(Data2, col=labels, main='Data type 2')
}
