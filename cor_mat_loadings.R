# Code by Ib Thorsgaard Jensen

# Loadings method for construction of correlation matrix
random_cor_mat <- function(d, k){
  W <- matrix(rnorm(d*k), d, k)
  S <- W %*% t(W) + diag(abs(rnorm(d)))
  C <- cov2cor(S)
  l <- which.max(rowSums(abs(C)))
  r <- 1:d; r[l] <- d; r[d] <- l
  C <- C[r,r]
  return(C)
}
