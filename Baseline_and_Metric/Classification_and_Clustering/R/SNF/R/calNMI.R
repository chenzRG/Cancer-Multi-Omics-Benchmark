calNMI <- function(x, y) {
##This function calculate the NMI between two clusters.	
	
	x = as.vector(x);
	y = as.vector(y);  
    return(max(0, .mutualInformation(x, y)/sqrt(.entropy(x) * .entropy(y)), na.rm=TRUE))
}


