\name{standardNormalization}
\alias{standardNormalization}
\title{
Standard Normalization
}
\description{
Normalize each column of the input data to have mean 0 and standard deviation 1.
}
\usage{
standardNormalization(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The unnormalized data.
}
}
\value{
The data normalized. 
}
\author{
Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

## Data1 is of size n x d_1, 
## where n is the number of patients, d_1 is the number of genes, 
## Data2 is of size n x d_2, 
## where n is the number of patients, d_2 is the number of methylation
data(Data1)
data(Data2)

Data1 = standardNormalization(Data1);
Data2 = standardNormalization(Data2);
}
