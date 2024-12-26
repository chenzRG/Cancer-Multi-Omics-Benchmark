chiDist2 <- function(A,B){
	###This function implements the Chi-Square distance between A and B
  n = nrow(A)
  m = nrow(B)
   d = ncol(A)
   stopifnot(d == ncol(B))
   res = matrix(nrow = n, ncol = m)
   sqA = A^2
   sqB = B^2
   twoAB = 2 * (A %*% t(B))
   for (a_num in 1:n) {
   	for (b_num in 1:m) {
   		res[a_num, b_num] = sum((sqA[a_num, ] + sqB[b_num, ] - twoAB[a_num, b_num]) / (A[a_num, ] + B[b_num, ]))
   		}
   		}
   		res = res / 2
   		return(res)  
}
