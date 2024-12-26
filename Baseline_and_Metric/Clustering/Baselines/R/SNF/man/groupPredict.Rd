\name{groupPredict}
\alias{groupPredict}
\title{
Group Predict
}
\description{
This function is used to predict the subtype of new patients.
}
\usage{
groupPredict(train, test, groups, K=20, alpha=0.5, t=20, method=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train}{
Training data. Has the same number of view and columns as test data. 
}
  \item{test}{
Test data. Has the same number of view and columns as training data. 
}
  \item{groups}{
The label for the training data.
}
  \item{K}{
Number of neighbors.
}
  \item{alpha}{
Hyperparameter used in constructing similarity network.
}
  \item{t}{
Number of iterations.
}
  \item{method}{
A indicator of which method to use to predict the label. method = 0 means to use local and global consistency; method = 1 means to use label propagation.
}
}
\value{
Returns the prediction of which group the test data belongs to. 
}
\author{
Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{

# Provide an example of predicting the new labels with label propagation

# Load views into list "dataL" and the cluster assignment into vector "label"
data(dataL)
data(label)

# Create the training and test data
n = floor(0.8*length(label)) # number of training cases
trainSample = sample.int(length(label), n)
train = lapply(dataL, function(x) x[trainSample, ]) # Use the first 150 samples for training
test = lapply(dataL, function(x) x[-trainSample, ]) # Test the rest of the data set
groups = label[trainSample]

# Set the other
K = 20
alpha = 0.5
t = 20
method = TRUE

# Apply the prediction function to the data
newLabel = groupPredict(train,test,groups,K,alpha,t,method)

# Compare the prediction accuracy
accuracy = sum(label[-trainSample] == newLabel[-c(1:n)])/(length(label) - n)

}
