displayClusters <- function(W, group) {
  normalize <- function(X) X / rowSums(X)
  ind = sort(as.vector(group),index.return=TRUE)
  ind = ind$ix
  diag(W) = 0
  W = normalize(W);
  W = W + t(W);
  image(1:ncol(W),1:nrow(W),W[ind,ind],col = grey(100:0 / 100),xlab = 'Patients',ylab='Patients');
}
